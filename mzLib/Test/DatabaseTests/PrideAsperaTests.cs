using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
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

/// <summary>
/// W2: one Aspera attempt (<see cref="PrideArchiveClient.DownloadFileViaAsperaAsync"/>), the opt-in policy over Aspera
/// and HTTPS (<see cref="PrideArchiveClient.DownloadFileWithPolicyAsync"/>), and the process-wide Aspera route. ascp is
/// replaced by a scripted session, so none of this needs ascp installed.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideAsperaTests
{
    // ---- test doubles -------------------------------------------------------

    /// <summary>A scripted ascp: it exits when told to, and records being stopped.</summary>
    private sealed class FakeAscp : IAscpSession
    {
        private readonly TaskCompletionSource<int> _exit = new(TaskCreationOptions.RunContinuationsAsynchronously);
        public Task<int> Exited => _exit.Task;
        public string LastLine { get; private set; } = string.Empty;
        public int Kills { get; private set; }

        public void Kill()
        {
            Kills++;
            _exit.TrySetResult(-1);
        }

        public FakeAscp Exit(int code, string line = "")
        {
            LastLine = line;
            _exit.TrySetResult(code);
            return this;
        }

        public void Dispose() { }
    }

    private sealed class StubHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, HttpResponseMessage> _responder;
        public List<HttpRequestMessage> Requests { get; } = new();

        public StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) => _responder = responder;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests.Add(request);
            return Task.FromResult(_responder(request));
        }
    }

    // ---- fixtures -------------------------------------------------------------

    private const string AsperaPath = "pride/data/archive/2018/10/PXD010820/run1.raw";
    private const string AsperaLocation = "prd_ascp@fasp.ebi.ac.uk:" + AsperaPath; // as PRIDE writes it: no '/' after ':'
    private const string HttpsUrl = "https://ftp.pride.ebi.ac.uk/" + AsperaPath;
    private const string FtpUrl = "ftp://ftp.pride.ebi.ac.uk/" + AsperaPath;

    /// <summary>The SHA-256 of EBI's public Aspera key, as Resources/PROVENANCE.md records it.</summary>
    private const string KeySha256 = "4f90eaab1230f9793d357ea7669ff9205edd30864a0bcf5ef3eb42a2720de1f1";

    private static readonly byte[] Body = Encoding.ASCII.GetBytes("ten bytes!");
    private static string Md5Of(byte[] bytes) => Convert.ToHexString(MD5.HashData(bytes)).ToLowerInvariant();
    private static readonly PrideFileChecksum Expected = new("run1.raw", Md5Of(Body), Body.Length);

    private static PrideArchiveFile MakeFile(params (string accession, string value)[] locations) =>
        new()
        {
            FileName = "run1.raw",
            FileSizeBytes = 999_999,
            FileCategory = new CvParam("PRIDE", "PRIDE:0000404", "category", "RAW"),
            PublicFileLocations = locations.Select(l => new CvParam("PRIDE", l.accession, "location", l.value)).ToList()
        };

    private static readonly (string, string) Aspera = (PrideArchiveExtensions.AsperaLocationAccession, AsperaLocation);
    private static readonly (string, string) Ftp = (PrideArchiveExtensions.FtpLocationAccession, FtpUrl);

    private static readonly PrideAsperaOptions Options = new() { AscpPath = @"C:\nowhere\ascp.exe", StallTimeout = TimeSpan.FromSeconds(30) };

    private string _tempDir;
    private string Destination => Path.Combine(_tempDir, "run1.raw");
    private string AsperaFolder => Path.Combine(_tempDir, PrideArchiveClient.AsperaFolderName);
    private List<ProcessStartInfo> _started;
    private List<FakeAscp> _sessions;
    private List<TimeSpan> _delays;
    private DateTimeOffset _now;

    [SetUp]
    public void SetUp()
    {
        _tempDir = Path.Combine(Path.GetTempPath(), "PrideAsperaTests", Guid.NewGuid().ToString("N"));
        _started = new();
        _sessions = new();
        _delays = new();
        _now = new DateTimeOffset(2026, 11, 16, 12, 0, 0, TimeSpan.Zero);
    }

    [TearDown]
    public void TearDown()
    {
        try { if (Directory.Exists(_tempDir)) Directory.Delete(_tempDir, recursive: true); }
        catch (IOException) { /* best effort */ }
    }

    /// <summary>
    /// A client whose ascp is <paramref name="ascp"/>, whose waits are recorded (and move the fake clock) instead of
    /// slept, and whose Aspera route is its own, on the fake clock.
    /// </summary>
    private PrideArchiveClient ClientWith(Func<ProcessStartInfo, FakeAscp> ascp, StubHandler handler = null,
        AsperaRoute route = null, int maxRetries = 0) =>
        new(new HttpClient(handler ?? new StubHandler(_ => throw new AssertionException("no HTTPS request was expected"))))
        {
            MaxRetries = maxRetries,
            RetryDelay = (wait, _) => { _delays.Add(wait); _now += wait; return Task.CompletedTask; },
            StartAscp = psi =>
            {
                _started.Add(psi);
                FakeAscp session = ascp(psi);
                _sessions.Add(session);
                return session;
            },
            AsperaRouteState = route ?? new AsperaRoute(() => _now)
        };

    private static string StagingOf(ProcessStartInfo psi) => psi.ArgumentList[^1];

    /// <summary>ascp writes <paramref name="bytes"/> under the source's name and exits with <paramref name="exit"/>.</summary>
    private static Func<ProcessStartInfo, FakeAscp> Writes(byte[] bytes, int exit = 0, string line = "Completed: 1 file") => psi =>
    {
        File.WriteAllBytes(Path.Combine(StagingOf(psi), "run1.raw"), bytes);
        return new FakeAscp().Exit(exit, line);
    };

    /// <summary>ascp writes nothing and exits with <paramref name="exit"/>, having said <paramref name="line"/> last.</summary>
    private static Func<ProcessStartInfo, FakeAscp> Says(int exit, string line) => _ => new FakeAscp().Exit(exit, line);

    private static HttpResponseMessage Download(byte[] body)
    {
        var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = new ByteArrayContent(body) };
        response.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
        return response;
    }

    // ---- TryGetAsperaSource ---------------------------------------------------

    [Test]
    public void TryGetAsperaSource_PridesForm_GainsTheSlashEbiDocuments()
    {
        Assert.That(MakeFile(Ftp, Aspera).TryGetAsperaSource(out string source), Is.True);
        Assert.That(source, Is.EqualTo("prd_ascp@fasp.ebi.ac.uk:/" + AsperaPath));
    }

    [Test]
    public void TryGetAsperaSource_WithTheSlashAlready_KeepsOne()
    {
        var file = MakeFile((PrideArchiveExtensions.AsperaLocationAccession, "prd_ascp@FASP.ebi.ac.uk:/" + AsperaPath));
        Assert.That(file.TryGetAsperaSource(out string source), Is.True);
        Assert.That(source, Is.EqualTo("prd_ascp@fasp.ebi.ac.uk:/" + AsperaPath));
    }

    [TestCase("prd_ascp@fasp.example.org:pride/run1.raw", Description = "another host: the bundled key opens EBI's only")]
    [TestCase("-oProxy=x@fasp.ebi.ac.uk:pride/run1.raw", Description = "a user ascp would read as an option")]
    [TestCase("prd ascp@fasp.ebi.ac.uk:pride/run1.raw", Description = "a user that is not a plain name")]
    [TestCase("prd_ascp@fasp.ebi.ac.uk:pride/../../etc/run1.raw", Description = "a '..' segment")]
    [TestCase("prd_ascp@fasp.ebi.ac.uk:pride\\run1.raw", Description = "a backslash")]
    [TestCase("prd_ascp@fasp.ebi.ac.uk:", Description = "no path")]
    [TestCase("prd_ascp@fasp.ebi.ac.uk:pride/", Description = "a folder, not a file")]
    [TestCase("fasp.ebi.ac.uk:pride/run1.raw", Description = "no user")]
    [TestCase("prd_ascp@fasp.ebi.ac.uk", Description = "no ':'")]
    public void TryGetAsperaSource_UnusableValue_IsSkipped(string value)
    {
        var file = MakeFile((PrideArchiveExtensions.AsperaLocationAccession, value));
        Assert.That(file.TryGetAsperaSource(out string source), Is.False);
        Assert.That(source, Is.Null);
    }

    [Test]
    public void TryGetAsperaSource_SelectsByAccession_NotByLookingLikeAspera()
    {
        Assert.That(MakeFile((PrideArchiveExtensions.FtpLocationAccession, AsperaLocation)).TryGetAsperaSource(out _), Is.False);
        Assert.That(MakeFile(Ftp).TryGetAsperaSource(out _), Is.False);
        Assert.That(new PrideArchiveFile { FileName = "run1.raw" }.TryGetAsperaSource(out _), Is.False);
        Assert.Throws<ArgumentNullException>(() => ((PrideArchiveFile)null).TryGetAsperaSource(out _));
    }

    // ---- the bundled key --------------------------------------------------------

    [Test]
    public void BundledKey_IsThePinnedEbiKey()
    {
        using Stream key = typeof(PrideArchiveClient).Assembly.GetManifestResourceStream(PrideArchiveClient.AsperaKeyResourceName);
        Assert.That(key, Is.Not.Null, "the key must be embedded in UsefulProteomicsDatabases.dll");
        Assert.That(Convert.ToHexString(SHA256.HashData(key)).ToLowerInvariant(), Is.EqualTo(KeySha256),
            "the embedded bytes must be the ones Resources/PROVENANCE.md pins (a CRLF checkout would change them)");
    }

    // ---- DownloadFileViaAsperaAsync: success ---------------------------------

    [Test]
    public async Task ViaAspera_Success_RunsEbisCommand_ChecksAndMovesTheCopy_AndLeavesNoFolder()
    {
        string keyDuringRun = null;
        using var client = ClientWith(psi =>
        {
            keyDuringRun = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(psi.ArgumentList[6]))).ToLowerInvariant();
            return Writes(Body)(psi);
        });

        string path = await client.DownloadFileViaAsperaAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, Options);

        Assert.That(path, Is.EqualTo(Destination));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
        ProcessStartInfo psi = _started.Single();
        Assert.That(psi.FileName, Is.EqualTo(Options.AscpPath));
        Assert.That(psi.UseShellExecute, Is.False);
        Assert.That(psi.ArgumentList.Take(6), Is.EqualTo(new[] { "-Q", "-l", "300m", "-P", "33001", "-i" }));
        Assert.That(psi.ArgumentList[7], Is.EqualTo("prd_ascp@fasp.ebi.ac.uk:/" + AsperaPath));
        Assert.That(StagingOf(psi), Does.StartWith(AsperaFolder + Path.DirectorySeparatorChar));
        Assert.That(psi.ArgumentList, Has.None.EqualTo("-T"), "encryption stays on: -T gave a damaged file");
        Assert.That(keyDuringRun, Is.EqualTo(KeySha256), "ascp reads the bundled key from a file while it runs");
        Assert.That(Directory.Exists(AsperaFolder), Is.False, "each attempt's folder, key included, is deleted");
    }

    [Test]
    public async Task ViaAspera_KeyPath_IsPassedAsGiven()
    {
        Directory.CreateDirectory(_tempDir);
        string myKey = Path.Combine(_tempDir, "my_key");
        File.WriteAllText(myKey, "key");
        using var client = ClientWith(Writes(Body));

        await client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected,
            new PrideAsperaOptions { AscpPath = Options.AscpPath, KeyPath = myKey, MaxRate = "50m" });

        Assert.That(_started.Single().ArgumentList[2], Is.EqualTo("50m"));
        Assert.That(_started.Single().ArgumentList[6], Is.EqualTo(myKey));
        Assert.That(File.Exists(myKey), Is.True, "a caller's key is never deleted");
    }

    [Test]
    public async Task ViaAspera_Success_RemovesAnOutOfDatePartialFromHttps()
    {
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination + ".partial", new byte[3]);
        File.WriteAllText(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix, "\"v1\"");
        using var client = ClientWith(Writes(Body));

        await client.DownloadFileViaAsperaAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, Options);

        Assert.That(File.Exists(Destination + ".partial"), Is.False);
        Assert.That(File.Exists(Destination + ".partial" + PrideArchiveClient.ValidatorSuffix), Is.False);
    }

    [Test]
    public async Task ViaAspera_NoOverwrite_GoodFileAlreadyThere_DoesNotRunAscp()
    {
        Directory.CreateDirectory(_tempDir);
        File.WriteAllBytes(Destination, Body);
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"));

        string path = await client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options, overwrite: false, verifyMd5: true);

        Assert.That(path, Is.EqualTo(Destination));
        Assert.That(_started, Is.Empty);
    }

    // ---- DownloadFileViaAsperaAsync: the copy's check ----------------------------

    [Test]
    public void ViaAspera_WrongSize_IsADamagedCopy_DiscardedAndTransient()
    {
        using var client = ClientWith(Writes(Encoding.ASCII.GetBytes("eleven byte")));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.DamagedCopy));
        Assert.That(e.AscpExitCode, Is.EqualTo(0));
        Assert.That(e.StatusCode, Is.Null);
        Assert.That(ServiceFailure.IsTransient(e), Is.True);
        Assert.That(e.Message, Does.Contain("'run1.raw' from fasp.ebi.ac.uk came to 11 bytes"));
        Assert.That(e.InnerException, Is.InstanceOf<InvalidDataException>());
        Assert.That(File.Exists(Destination), Is.False);
        Assert.That(Directory.Exists(AsperaFolder), Is.False);
    }

    [Test]
    public void ViaAspera_EbisErrorText_WithoutMd5_IsADamagedCopy()
    {
        byte[] damaged = Encoding.ASCII.GetBytes("<Error><Code>X");
        using var client = ClientWith(Writes(damaged));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() => client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(Body), damaged.Length), Options));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.DamagedCopy));
        Assert.That(e.Message, Does.Contain("<Error><Code>"));
    }

    [Test]
    public void ViaAspera_EbisExactOpening_IsRefusedEvenWhenTheMd5Matches()
    {
        // Promised in PXReprise 033, and what an HTTPS download stopped at the opening does.
        byte[] damaged = Encoding.ASCII.GetBytes("raw<Error><Code>ConnectionClosedException</Code>raw");
        using var client = ClientWith(Writes(damaged));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() => client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir,
            new PrideFileChecksum("run1.raw", Md5Of(damaged), damaged.Length), Options, verifyMd5: true));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.DamagedCopy));
        Assert.That(e.Message, Does.Contain("at byte 3"));
        Assert.That(File.Exists(Destination), Is.False);
    }

    [Test]
    public void ViaAspera_WrongMd5_IsADamagedCopy()
    {
        byte[] other = Encoding.ASCII.GetBytes("TEN BYTES!");
        using var client = ClientWith(Writes(other));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options, verifyMd5: true));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.DamagedCopy));
        Assert.That(e.Message, Does.Contain("has MD5 " + Md5Of(other)));
    }

    // ---- DownloadFileViaAsperaAsync: ascp's failures ----------------------------

    // The lines are ascp's own, as PXReprise's batch logged them, 2026-10-08 to 2026-11-15 (counts in brackets).
    [TestCase(1, "Session Stop  (Error: failed to authenticate)", PrideAsperaFailure.LoginRefused)] // 2,272
    [TestCase(1, "Session Stop  (Error: Unable to connect via SSH)", PrideAsperaFailure.CannotConnect)] // 468
    [TestCase(1, "Session Stop  (Error: Client unable to connect to server (check UDP port and firewall))", PrideAsperaFailure.CannotConnect)] // 12
    [TestCase(1, "Session Stop  (Error: Failed to open TCP connection for SSH)", PrideAsperaFailure.CannotConnect)] // 1
    [TestCase(1, "Session Stop  (Error: Session initiation failed, read timed out)", PrideAsperaFailure.CannotConnect)] // 1
    [TestCase(1, "Session Stop  (Error: Session data transfer timeout)", PrideAsperaFailure.AscpFailed)] // 91
    [TestCase(1, "Session Stop  (Error: Session shutdown failed, Failure while draining incoming flow)", PrideAsperaFailure.AscpFailed)] // 13
    [TestCase(1, "Session Stop  (Error: Server aborted session: Failed to connect to management socket)", PrideAsperaFailure.AscpFailed)] // 3
    [TestCase(-1, "", PrideAsperaFailure.AscpFailed)] // 1
    public void ViaAspera_AscpExits_FailureIsClassifiedFromWhatItSaid(int exit, string line, PrideAsperaFailure failure)
    {
        using var client = ClientWith(Says(exit, line));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options));

        Assert.That(e.Failure, Is.EqualTo(failure));
        Assert.That(e.AscpExitCode, Is.EqualTo(exit));
        Assert.That(e.StatusCode, Is.Null);
        Assert.That(ServiceFailure.IsTransient(e), Is.True, "every Aspera failure is worth trying again later");
        Assert.That(e.InnerException, Is.InstanceOf<IOException>(), "PXReprise's own check needs an inner exception");
        Assert.That(e.Message, Does.StartWith("The PRIDE Aspera download of 'run1.raw' from fasp.ebi.ac.uk"));
        Assert.That(e.Message, Does.Contain($"ascp exited with {exit}"));
        Assert.That(Directory.Exists(AsperaFolder), Is.False);
    }

    [Test]
    public void ViaAspera_DiskFull_IsTheCallersDisk_NotTransient()
    {
        using var client = ClientWith(Says(1, "Session Stop  (Error: No space left on disk)")); // 54 in the batch

        var e = Assert.ThrowsAsync<IOException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options));

        Assert.That(e, Is.Not.InstanceOf<PrideAsperaException>());
        Assert.That(ServiceFailure.IsTransient(e), Is.False);
        Assert.That(e.Message, Does.Contain("ran out of disk space"));
    }

    [Test]
    public void ViaAspera_ExitsZeroWithNoFile_IsAscpFailed()
    {
        using var client = ClientWith(Says(0, "Completed: 0 files"));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.AscpFailed));
        Assert.That(e.AscpExitCode, Is.EqualTo(0));
    }

    [Test]
    public void ViaAspera_MessagesNameNeitherTheSourcePathNorLocalFolders()
    {
        using var client = ClientWith(psi => new FakeAscp().Exit(1,
            $"Session Stop  (Error: failed for {psi.ArgumentList[7]} into {StagingOf(psi)} with {psi.ArgumentList[6]})"));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options));

        string everything = e.ToString();
        Assert.That(everything, Does.Not.Contain("pride/data/archive"));
        Assert.That(everything, Does.Not.Contain(_tempDir));
        Assert.That(everything, Does.Contain("<source>").And.Contain("<folder>"));
    }

    [Test]
    public void ViaAspera_NoBytesForTheStallTimeout_StopsAscp()
    {
        using var client = ClientWith(_ => new FakeAscp()); // never exits, never writes

        var stopwatch = Stopwatch.StartNew();
        var e = Assert.ThrowsAsync<PrideAsperaException>(() => client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir,
            Expected, new PrideAsperaOptions { AscpPath = Options.AscpPath, StallTimeout = TimeSpan.FromMilliseconds(300) }));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.Stalled));
        Assert.That(e.AscpExitCode, Is.Null);
        Assert.That(e.InnerException, Is.InstanceOf<TimeoutException>());
        Assert.That(ServiceFailure.IsTransient(e), Is.True);
        Assert.That(_sessions.Single().Kills, Is.EqualTo(1));
        Assert.That(stopwatch.Elapsed, Is.LessThan(TimeSpan.FromSeconds(10)));
        Assert.That(Directory.Exists(AsperaFolder), Is.False);
    }

    [Test]
    public async Task ViaAspera_SlowButReceiving_IsNeverStopped()
    {
        // Bytes arrive every 50 ms for 1.2 s, four times the 300 ms stall window: only silence counts.
        using var client = ClientWith(psi =>
        {
            var session = new FakeAscp();
            string inProgress = Path.Combine(StagingOf(psi), "run1.raw.aspx");
            _ = Task.Run(async () =>
            {
                for (int i = 1; i <= 24; i++)
                {
                    await Task.Delay(50);
                    File.WriteAllBytes(inProgress, new byte[i]);
                }
                File.Delete(inProgress);
                File.WriteAllBytes(Path.Combine(StagingOf(psi), "run1.raw"), Body);
                session.Exit(0);
            });
            return session;
        });

        await client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected,
            new PrideAsperaOptions { AscpPath = Options.AscpPath, StallTimeout = TimeSpan.FromMilliseconds(300) });

        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
        Assert.That(_sessions.Single().Kills, Is.Zero);
    }

    [Test]
    public void ViaAspera_Cancelled_StopsAscpAndLeavesNoFolder()
    {
        using var client = ClientWith(_ => new FakeAscp());
        using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(100));

        Assert.CatchAsync<OperationCanceledException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options, cancellationToken: cts.Token));

        Assert.That(_sessions.Single().Kills, Is.EqualTo(1));
        Assert.That(Directory.Exists(AsperaFolder), Is.False);
    }

    [Test]
    public void ViaAspera_AscpCannotStart_IsAConfigurationError_NotAnOutage()
    {
        using var client = ClientWith(_ => throw new System.ComponentModel.Win32Exception(2, "The system cannot find the file specified."));

        var e = Assert.ThrowsAsync<MzLibException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected, Options));

        Assert.That(e.Message, Does.Contain(Options.AscpPath));
        Assert.That(ServiceFailure.IsTransient(e), Is.False);
        Assert.That(Directory.Exists(AsperaFolder), Is.False);
    }

    [Test]
    public void ViaAspera_RealStart_OfAMissingExecutable_IsAConfigurationError()
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => throw new AssertionException("no HTTPS"))));

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected,
            new PrideAsperaOptions { AscpPath = Path.Combine(_tempDir, "no-such-ascp.exe") }));

        Assert.That(e.Message, Does.StartWith("ascp could not be started"));
    }

    [Test]
    public void ViaAspera_MissingKeyPath_IsAConfigurationError()
    {
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"));

        Assert.ThrowsAsync<MzLibException>(() => client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected,
            new PrideAsperaOptions { AscpPath = Options.AscpPath, KeyPath = Path.Combine(_tempDir, "missing_key") }));
    }

    [Test]
    public void ViaAspera_NoAsperaLocation_IsNotSupported()
    {
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"));

        var e = Assert.ThrowsAsync<NotSupportedException>(() =>
            client.DownloadFileViaAsperaAsync(MakeFile(Ftp), _tempDir, Expected, Options));
        Assert.That(e.Message, Does.Contain("no Aspera location"));
    }

    [Test]
    public void ViaAspera_BadArguments_AreRefusedBeforeAscpRuns()
    {
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"));
        PrideArchiveFile file = MakeFile(Aspera);

        Assert.ThrowsAsync<ArgumentNullException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir, null, Options));
        Assert.ThrowsAsync<ArgumentNullException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir, Expected, null));
        Assert.ThrowsAsync<ArgumentException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir, Expected, new PrideAsperaOptions()));
        foreach (string rate in new[] { "-T", "300x", "m", "", " 300m" })
            Assert.ThrowsAsync<ArgumentException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir, Expected,
                new PrideAsperaOptions { AscpPath = Options.AscpPath, MaxRate = rate }), rate);
        Assert.ThrowsAsync<ArgumentOutOfRangeException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir, Expected,
            new PrideAsperaOptions { AscpPath = Options.AscpPath, StallTimeout = TimeSpan.Zero }));
        Assert.ThrowsAsync<ArgumentException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir,
            new PrideFileChecksum("run2.raw", Md5Of(Body), Body.Length), Options));
        Assert.That(_started, Is.Empty);
    }

    // ---- the real process wrapper -------------------------------------------------

    [Test]
    public async Task AscpProcess_ReportsTheExitCodeAndTheLastLine()
    {
        using IAscpSession session = AscpProcess.Start(new ProcessStartInfo("cmd.exe", "/c echo first 1>&2 & echo Session Stop (Error: last) 1>&2 & exit 3")
        {
            UseShellExecute = false, RedirectStandardOutput = true, RedirectStandardError = true, CreateNoWindow = true
        });

        Assert.That(await session.Exited.WaitAsync(TimeSpan.FromSeconds(30)), Is.EqualTo(3));
        Assert.That(session.LastLine, Is.EqualTo("Session Stop (Error: last)"));
    }

    [Test]
    public async Task AscpProcess_Kill_EndsTheProcessTree()
    {
        using IAscpSession session = AscpProcess.Start(new ProcessStartInfo("cmd.exe", "/c ping -n 60 127.0.0.1 > nul")
        {
            UseShellExecute = false, RedirectStandardOutput = true, RedirectStandardError = true, CreateNoWindow = true
        });

        session.Kill();

        await session.Exited.WaitAsync(TimeSpan.FromSeconds(30));
        session.Kill(); // a second stop does nothing
    }

    // ---- DownloadFileWithPolicyAsync -------------------------------------------------

    [Test]
    public void Policy_Defaults_ArePXReprisesValues()
    {
        var policy = new PrideDownloadPolicy();
        Assert.That(policy.Aspera, Is.Null);
        Assert.That(policy.Attempts, Is.EqualTo(8));
        Assert.That(policy.Backoff, Is.EqualTo(TimeSpan.FromSeconds(15)));
        Assert.That(policy.AsperaTries, Is.EqualTo(3));
        Assert.That(policy.HttpsAfterDamagedAsperaCopy, Is.True);
        Assert.That(policy.AsperaPauseAfterRouteFailure, Is.EqualTo(TimeSpan.FromMinutes(30)));
        Assert.That(policy.MaxConcurrentAsperaSessions, Is.EqualTo(3));
        Assert.That(policy.AsperaStartGap, Is.EqualTo(TimeSpan.FromSeconds(10)));
        Assert.That(new PrideAsperaOptions().MaxRate, Is.EqualTo("300m"));
        Assert.That(new PrideAsperaOptions().StallTimeout, Is.EqualTo(TimeSpan.FromMinutes(5)));
    }

    [Test]
    public async Task Policy_AsperaTriesThenHttps_WithGrowingWaits_ReportingEachTry()
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientWith(Says(1, "Session Stop  (Error: Session data transfer timeout)"), handler);
        var tries = new List<PrideDownloadTry>();

        string path = await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected,
            new PrideDownloadPolicy { Aspera = Options, Tries = new SyncProgress<PrideDownloadTry>(tries.Add) });

        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Body));
        Assert.That(_started, Has.Count.EqualTo(3));
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
        Assert.That(_delays, Is.EqualTo(new[] { 15, 30, 45 }.Select(s => TimeSpan.FromSeconds(s))));
        Assert.That(tries.Select(t => t.Route), Is.EqualTo(new[]
            { PrideDownloadRoute.Aspera, PrideDownloadRoute.Aspera, PrideDownloadRoute.Aspera, PrideDownloadRoute.Https }));
        Assert.That(tries.Select(t => t.Outcome), Is.EqualTo(new[]
            { PrideDownloadOutcome.Failed, PrideDownloadOutcome.Failed, PrideDownloadOutcome.Failed, PrideDownloadOutcome.Accepted }));
        Assert.That(tries.Select(t => t.Try), Is.EqualTo(new[] { 1, 2, 3, 4 }));
        Assert.That(tries.Select(t => t.Host), Is.EqualTo(new[] { "fasp.ebi.ac.uk", "fasp.ebi.ac.uk", "fasp.ebi.ac.uk", "ftp.pride.ebi.ac.uk" }));
        Assert.That(tries[0].AsperaFailure, Is.EqualTo(PrideAsperaFailure.AscpFailed));
        Assert.That(tries[3].Error, Is.Null);
        Assert.That(tries.All(t => t.FileName == "run1.raw" && t.FileSizeBytes == Body.Length), Is.True);
    }

    [Test]
    public async Task Policy_LoginRefused_SetsAsperaAsideForEveryFile_ThenTestsItAgain()
    {
        var route = new AsperaRoute(() => _now);
        var handler = new StubHandler(_ => Download(Body));
        using var refused = ClientWith(Says(1, "Session Stop  (Error: failed to authenticate)"), handler, route);
        var policy = new PrideDownloadPolicy { Aspera = Options };

        await refused.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, policy);
        Assert.That(_started, Has.Count.EqualTo(1), "one refusal; the next try went over HTTPS");
        Assert.That(route.IsPaused, Is.True);

        // Another client in the same process: the route is still set aside, so it goes straight to HTTPS.
        File.Delete(Destination);
        using var other = ClientWith(Writes(Body), handler, route);
        await other.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, policy);
        Assert.That(_started, Has.Count.EqualTo(1));

        // Thirty minutes on, the first try tests Aspera again.
        _now += TimeSpan.FromMinutes(30);
        File.Delete(Destination);
        await other.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, policy);
        Assert.That(_started, Has.Count.EqualTo(2));
        Assert.That(route.IsPaused, Is.False);
    }

    [Test]
    public async Task Policy_DamagedAsperaCopy_SendsTheRestToHttps()
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientWith(Writes(Encoding.ASCII.GetBytes("eleven byte")), handler);
        var tries = new List<PrideDownloadTry>();

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected,
            new PrideDownloadPolicy { Aspera = Options, Tries = new SyncProgress<PrideDownloadTry>(tries.Add) });

        Assert.That(_started, Has.Count.EqualTo(1));
        Assert.That(tries.Select(t => t.Outcome), Is.EqualTo(new[] { PrideDownloadOutcome.Damaged, PrideDownloadOutcome.Accepted }));
        Assert.That(tries[1].Route, Is.EqualTo(PrideDownloadRoute.Https));
    }

    [Test]
    public async Task Policy_DamagedAsperaCopy_StaysOnAspera_WhenAskedTo()
    {
        int run = 0;
        using var client = ClientWith(psi => Writes(++run == 1 ? Encoding.ASCII.GetBytes("eleven byte") : Body)(psi),
            new StubHandler(_ => throw new AssertionException("no HTTPS request was expected")));

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected,
            new PrideDownloadPolicy { Aspera = Options, HttpsAfterDamagedAsperaCopy = false });

        Assert.That(_started, Has.Count.EqualTo(2));
    }

    [Test]
    public void Policy_AsperaOnlyFile_KeepsToAspera_EvenWhileItIsSetAside()
    {
        using var client = ClientWith(Says(1, "Session Stop  (Error: failed to authenticate)"));

        var e = Assert.ThrowsAsync<PrideAsperaException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Aspera), _tempDir,
            Expected, new PrideDownloadPolicy { Aspera = Options, Attempts = 4 }));

        Assert.That(e.Failure, Is.EqualTo(PrideAsperaFailure.LoginRefused));
        Assert.That(_started, Has.Count.EqualTo(4), "with no HTTPS location, Aspera is the only route");
    }

    [Test]
    public void Policy_HttpsTry_IsOneAttempt_SoTheBudgetsDoNotMultiply()
    {
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.ServiceUnavailable));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler, maxRetries: 3);

        var e = Assert.ThrowsAsync<HttpRequestException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir,
            Expected, new PrideDownloadPolicy { Attempts = 2 }));

        Assert.That(e.StatusCode, Is.EqualTo(HttpStatusCode.ServiceUnavailable));
        Assert.That(ServiceFailure.IsTransient(e), Is.True);
        Assert.That(handler.Requests, Has.Count.EqualTo(2), "two tries, not two times the HTTPS download's four attempts");
        Assert.That(_delays, Is.EqualTo(new[] { TimeSpan.FromSeconds(15) }));
    }

    [Test]
    public void Policy_Refusal_IsThrownAtOnce()
    {
        var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.NotFound));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);
        var tries = new List<PrideDownloadTry>();

        var e = Assert.ThrowsAsync<HttpRequestException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir,
            Expected, new PrideDownloadPolicy { Tries = new SyncProgress<PrideDownloadTry>(tries.Add) }));

        Assert.That(e.StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
        Assert.That(handler.Requests, Has.Count.EqualTo(1));
        Assert.That(tries.Single().Outcome, Is.EqualTo(PrideDownloadOutcome.Failed));
    }

    [Test]
    public void Policy_DamagedHttpsCopies_AreTriedAgain_ThenReportedAsData()
    {
        var handler = new StubHandler(_ => Download(Encoding.ASCII.GetBytes("eleven byte")));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);
        var tries = new List<PrideDownloadTry>();

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir,
            Expected, new PrideDownloadPolicy { Attempts = 3, Tries = new SyncProgress<PrideDownloadTry>(tries.Add) }));

        Assert.That(e.Message, Does.Contain("came to 11 bytes"));
        Assert.That(handler.Requests, Has.Count.EqualTo(3));
        Assert.That(tries.Select(t => t.Outcome), Is.All.EqualTo(PrideDownloadOutcome.Damaged));
        Assert.That(File.Exists(Destination), Is.False);
    }

    [Test]
    public void Policy_AscpCannotStart_IsThrownAtOnce()
    {
        using var client = ClientWith(_ => throw new System.ComponentModel.Win32Exception(2, "missing"), new StubHandler(_ => Download(Body)));

        Assert.ThrowsAsync<MzLibException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir,
            Expected, new PrideDownloadPolicy { Aspera = Options }));
    }

    [Test]
    public async Task Policy_WithoutAspera_IsHttpsOnly_AndSkipsAGoodFileAlreadyThere()
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, new PrideDownloadPolicy());
        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, new PrideDownloadPolicy(), overwrite: false);

        Assert.That(handler.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public void Policy_NoRouteAtAll_IsNotSupported()
    {
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"));

        Assert.ThrowsAsync<NotSupportedException>(() =>
            client.DownloadFileWithPolicyAsync(MakeFile(Aspera), _tempDir, Expected, new PrideDownloadPolicy()));
        Assert.ThrowsAsync<NotSupportedException>(() =>
            client.DownloadFileWithPolicyAsync(MakeFile(), _tempDir, Expected, new PrideDownloadPolicy { Aspera = Options }));
    }

    [Test]
    public void Policy_BadSettings_AreRefused()
    {
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"));
        PrideArchiveFile file = MakeFile(Ftp);

        Assert.ThrowsAsync<ArgumentNullException>(() => client.DownloadFileWithPolicyAsync(file, _tempDir, Expected, null));
        foreach (var policy in new[]
        {
            new PrideDownloadPolicy { Attempts = 0 },
            new PrideDownloadPolicy { AsperaTries = -1 },
            new PrideDownloadPolicy { MaxConcurrentAsperaSessions = 0 },
            new PrideDownloadPolicy { Backoff = TimeSpan.FromSeconds(-1) },
            new PrideDownloadPolicy { AsperaStartGap = TimeSpan.FromSeconds(-1) },
            new PrideDownloadPolicy { AsperaPauseAfterRouteFailure = TimeSpan.FromSeconds(-1) },
            new PrideDownloadPolicy { Aspera = new PrideAsperaOptions() }
        })
            Assert.ThrowsAsync<ArgumentException>(() => client.DownloadFileWithPolicyAsync(file, _tempDir, Expected, policy));
    }

    // ---- AsperaRoute -------------------------------------------------------------------

    [Test]
    public async Task Route_Cap_HoldsTheNextSessionUntilOneEnds()
    {
        var route = new AsperaRoute(() => _now);
        Func<TimeSpan, CancellationToken, Task> delay = (wait, _) => { _now += wait; return Task.CompletedTask; };

        IDisposable first = await route.EnterAsync(1, TimeSpan.Zero, delay, CancellationToken.None);
        Task<IDisposable> second = route.EnterAsync(1, TimeSpan.Zero, delay, CancellationToken.None);
        await Task.Delay(50);
        Assert.That(second.IsCompleted, Is.False);

        first.Dispose();
        (await second.WaitAsync(TimeSpan.FromSeconds(10))).Dispose();
        first.Dispose(); // a second dispose does nothing

        IDisposable third = await route.EnterAsync(1, TimeSpan.Zero, delay, CancellationToken.None).WaitAsync(TimeSpan.FromSeconds(10));
        third.Dispose();
    }

    [Test]
    public async Task Route_StartGap_SpacesStarts()
    {
        var route = new AsperaRoute(() => _now);
        var waits = new List<TimeSpan>();
        Func<TimeSpan, CancellationToken, Task> delay = (wait, _) => { waits.Add(wait); _now += wait; return Task.CompletedTask; };

        using (await route.EnterAsync(3, TimeSpan.FromSeconds(10), delay, CancellationToken.None)) { }
        _now += TimeSpan.FromSeconds(4);
        using (await route.EnterAsync(3, TimeSpan.FromSeconds(10), delay, CancellationToken.None)) { }

        Assert.That(waits, Is.EqualTo(new[] { TimeSpan.FromSeconds(6) }));
    }

    [Test]
    public void Route_WaitingForASession_CanBeCancelled()
    {
        var route = new AsperaRoute(() => _now);
        using var cts = new CancellationTokenSource();
        Func<TimeSpan, CancellationToken, Task> delay = (_, _) => Task.CompletedTask;

        using IDisposable held = route.EnterAsync(1, TimeSpan.Zero, delay, CancellationToken.None).Result;
        Task<IDisposable> waiting = route.EnterAsync(1, TimeSpan.Zero, delay, cts.Token);
        cts.Cancel();

        Assert.CatchAsync<OperationCanceledException>(() => waiting);
    }

    [Test]
    public void Route_Pause_KeepsTheLongerOne_AndZeroSetsNothingAside()
    {
        var route = new AsperaRoute(() => _now);
        route.Pause(TimeSpan.Zero);
        Assert.That(route.IsPaused, Is.False);

        route.Pause(TimeSpan.FromMinutes(30));
        route.Pause(TimeSpan.FromMinutes(5));
        _now += TimeSpan.FromMinutes(29);
        Assert.That(route.IsPaused, Is.True);
        _now += TimeSpan.FromMinutes(1);
        Assert.That(route.IsPaused, Is.False);
    }

    /// <summary>An <see cref="IProgress{T}"/> that reports on the caller's thread, so a test can read the list at once.</summary>
    private sealed class SyncProgress<T>(Action<T> report) : IProgress<T>
    {
        public void Report(T value) => report(value);
    }
}

/// <summary>
/// Live canary for W2 against EBI's Aspera server: one small public file, checked against PRIDE's own checksum row.
/// Skips when no ascp is installed (CI has none; set <c>MZLIB_ASCP</c> to its path, or install IBM Aspera for Desktop),
/// and, through <see cref="ExternalServiceTestHelper.RunAsync"/>, when EBI fails: every Aspera failure is an
/// <see cref="HttpRequestException"/>. A configuration error or a contract break fails.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideAsperaLiveTests
{
    private const string DesktopAscp = @"C:\Program Files\IBM Aspera\transferd\bin\ascp.exe";

    [Test]
    public async Task DownloadFileViaAsperaAsync_LiveSmallFile_MatchesPridesChecksum()
    {
        string ascp = Environment.GetEnvironmentVariable("MZLIB_ASCP");
        if (string.IsNullOrWhiteSpace(ascp))
            ascp = DesktopAscp;
        if (!File.Exists(ascp))
            Assert.Ignore($"No ascp at '{ascp}'; set MZLIB_ASCP to run this canary.");

        await ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            // one attempt: an outage skips at once, not after the retry backoffs
            using var client = new PrideArchiveClient { MaxRetries = 0 };
            var files = await client.GetProjectFilesAsync("PXD000001");
            // Small, listed in both the manifest and the checksum list, and not an SDRF (PRIDE rewrites those in place).
            PrideArchiveFile mztab = files.Single(f => f.FileName == "F063721.dat-mztab.txt");
            var checksums = await client.GetFileChecksumsAsync("PXD000001");
            PrideFileChecksum expected = checksums["F063721.dat-mztab.txt"];

            string dir = Path.Combine(Path.GetTempPath(), "PrideAsperaLive", Guid.NewGuid().ToString("N"));
            try
            {
                string path = await client.DownloadFileViaAsperaAsync(mztab, dir, expected,
                    new PrideAsperaOptions { AscpPath = ascp }, verifyMd5: true);
                Assert.That(new FileInfo(path).Length, Is.EqualTo(expected.SizeBytes));
                Assert.That(Directory.Exists(Path.Combine(dir, PrideArchiveClient.AsperaFolderName)), Is.False);
            }
            finally
            {
                if (Directory.Exists(dir)) Directory.Delete(dir, recursive: true);
            }
        });
    }
}
