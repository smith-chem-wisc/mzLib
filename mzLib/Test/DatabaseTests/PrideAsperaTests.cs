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

    /// <summary>
    /// What a path argument names inside ascp's process. ascp starts in <see cref="ProcessStartInfo.WorkingDirectory"/>
    /// (which the OS takes from this process's directory, when relative), so a relative argument is relative to that
    /// folder, not to this process's. A rooted argument is unchanged.
    /// </summary>
    private static string InAscp(ProcessStartInfo psi, string argument) =>
        Path.Combine(Path.GetFullPath(psi.WorkingDirectory), argument);

    /// <summary>The folder ascp was told to write into, as ascp sees it.</summary>
    private static string StagingOf(ProcessStartInfo psi) => InAscp(psi, psi.ArgumentList[^1]);

    /// <summary>The key ascp was told to log in with (the argument after <c>-i</c>), as ascp sees it.</summary>
    private static string KeyOf(ProcessStartInfo psi) => InAscp(psi, psi.ArgumentList[psi.ArgumentList.IndexOf("-i") + 1]);

    /// <summary>
    /// ascp writes <paramref name="bytes"/> under the source's name and exits with <paramref name="exit"/>, as the real
    /// one would from its own folder: a key it cannot find there ends the session with ascp 4.4.8's own line (seen
    /// offline, 2026-11-19), and so does a destination folder it cannot find (a line of this fake's, since ascp
    /// reaches that check only after logging in).
    /// </summary>
    private static Func<ProcessStartInfo, FakeAscp> Writes(byte[] bytes, int exit = 0, string line = "Completed: 1 file") => psi =>
    {
        if (!File.Exists(KeyOf(psi)))
            return new FakeAscp().Exit(1, $"Session Stop  (Error: Private key file not found at path {KeyOf(psi)})");
        if (!Directory.Exists(StagingOf(psi)))
            return new FakeAscp().Exit(1, $"Session Stop  (Error: fake ascp: no folder at {StagingOf(psi)})");
        File.WriteAllBytes(Path.Combine(StagingOf(psi), "run1.raw"), bytes);
        return new FakeAscp().Exit(exit, line);
    };

    /// <summary>
    /// Points the test's folder at a new one under this process's current directory, and returns its relative name.
    /// (The temp folder can be on another drive than the current directory, where no relative path reaches it.)
    /// </summary>
    private string RelativeTempDir()
    {
        string relative = "PrideAsperaTests-" + Guid.NewGuid().ToString("N");
        _tempDir = Path.GetFullPath(relative);
        return relative;
    }

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
            keyDuringRun = Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(KeyOf(psi)))).ToLowerInvariant();
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

    // ascp starts in its staging folder (WorkingDirectory), so a relative path it is given names a file under that
    // folder, not under the caller's directory. These three are Alexander-Sol's cases on #1612.

    [Test]
    public async Task ViaAspera_RelativeDestination_IsResolvedBeforeAscpStarts_SoAscpFindsTheBundledKeyAndItsFolder()
    {
        string relative = RelativeTempDir();
        using var client = ClientWith(Writes(Body));

        string path = await client.DownloadFileViaAsperaAsync(MakeFile(Aspera), relative, Expected, Options);

        Assert.That(path, Is.EqualTo(Destination), "the full path, as documented");
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
        ProcessStartInfo psi = _started.Single();
        Assert.That(psi.WorkingDirectory, Does.StartWith(AsperaFolder + Path.DirectorySeparatorChar));
        Assert.That(psi.ArgumentList[6], Does.StartWith(AsperaFolder + Path.DirectorySeparatorChar), "the bundled key's full path");
        Assert.That(psi.ArgumentList[^1], Does.StartWith(AsperaFolder + Path.DirectorySeparatorChar), "the staging folder's full path");
        Assert.That(Directory.Exists(AsperaFolder), Is.False);
    }

    [Test]
    public async Task ViaAspera_RelativeDestination_WithAGivenKey_AscpWritesIntoItsFolder()
    {
        string relative = RelativeTempDir();
        Directory.CreateDirectory(_tempDir);
        string myKey = Path.Combine(_tempDir, "my_key");
        File.WriteAllText(myKey, "key");
        using var client = ClientWith(Writes(Body));

        string path = await client.DownloadFileViaAsperaAsync(MakeFile(Aspera), relative, Expected,
            new PrideAsperaOptions { AscpPath = Options.AscpPath, KeyPath = myKey });

        Assert.That(path, Is.EqualTo(Destination));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
        Assert.That(_started.Single().ArgumentList[^1], Does.StartWith(AsperaFolder + Path.DirectorySeparatorChar));
    }

    [Test]
    public async Task ViaAspera_RelativeKeyPath_WithAFullDestination_IsResolvedBeforeAscpStarts()
    {
        string relativeKey = Path.Combine(RelativeTempDir(), "my_key");
        Directory.CreateDirectory(_tempDir);
        File.WriteAllText(relativeKey, "key");
        using var client = ClientWith(Writes(Body));

        await client.DownloadFileViaAsperaAsync(MakeFile(Aspera), _tempDir, Expected,
            new PrideAsperaOptions { AscpPath = Options.AscpPath, KeyPath = relativeKey });

        Assert.That(_started.Single().ArgumentList[6], Is.EqualTo(Path.Combine(_tempDir, "my_key")));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
        Assert.That(File.Exists(relativeKey), Is.True, "a caller's key is never deleted");
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
        foreach (string key in new[] { "", " " })
            Assert.That(Assert.ThrowsAsync<ArgumentException>(() => client.DownloadFileViaAsperaAsync(file, _tempDir, Expected,
                new PrideAsperaOptions { AscpPath = Options.AscpPath, KeyPath = key })).Message,
                Does.Contain("PrideAsperaOptions.KeyPath"), $"KeyPath '{key}'");
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

    /// <summary>A shell running <paramref name="windows"/> or <paramref name="unix"/>, with ascp's redirections.</summary>
    private static ProcessStartInfo Shell(string windows, string unix)
    {
        var startInfo = OperatingSystem.IsWindows()
            ? new ProcessStartInfo("cmd.exe") { ArgumentList = { "/c", windows } }
            : new ProcessStartInfo("/bin/sh") { ArgumentList = { "-c", unix } };
        startInfo.UseShellExecute = false;
        startInfo.RedirectStandardOutput = true;
        startInfo.RedirectStandardError = true;
        startInfo.CreateNoWindow = true;
        return startInfo;
    }

    [Test]
    public async Task AscpProcess_ReportsTheExitCodeAndTheLastLine()
    {
        // Both lines on stderr: the two streams are read concurrently, so only one stream's order is defined.
        using IAscpSession session = AscpProcess.Start(Shell(
            "echo first 1>&2 & echo Session Stop (Error: last) 1>&2 & exit 3",
            "echo first 1>&2; echo 'Session Stop (Error: last)' 1>&2; exit 3"));

        Assert.That(await session.Exited.WaitAsync(TimeSpan.FromSeconds(30)), Is.EqualTo(3));
        Assert.That(session.LastLine, Is.EqualTo("Session Stop (Error: last)"));
    }

    [Test]
    public async Task AscpProcess_Kill_EndsTheProcessTree()
    {
        // The shell starts a child that would run for a minute and holds the redirected streams. Unless the child is
        // stopped too, Exited waits for it, past the 30 s below.
        using IAscpSession session = AscpProcess.Start(Shell("ping -n 60 127.0.0.1 > nul", "sleep 60; echo done"));

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

    // ---- PXReprise 040: the HTTPS location of each try, and what each try's record says ----------

    private const string GlobusHost = "globus.example.org";
    private const string GlobusUrl = "https://" + GlobusHost + "/" + AsperaPath;

    /// <summary>A body that delivers <paramref name="head"/>, then breaks as a dropped connection does.</summary>
    private sealed class BreaksAfter(byte[] head) : Stream
    {
        private int _at;

        public override ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            if (_at >= head.Length)
                return ValueTask.FromException<int>(new IOException("The response ended prematurely."));
            int n = Math.Min(buffer.Length, head.Length - _at);
            head.AsSpan(_at, n).CopyTo(buffer.Span);
            _at += n;
            return ValueTask.FromResult(n);
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

    /// <summary><see cref="Body"/>'s first <paramref name="received"/> bytes, then a dropped connection; resumable (ETag v1).</summary>
    private static HttpResponseMessage DropsAfter(int received)
    {
        var response = new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new BreaksAfter(Body[..received])) };
        response.Content.Headers.ContentLength = Body.Length;
        response.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
        return response;
    }

    /// <summary>The rest of <see cref="Body"/> from <paramref name="from"/>, as a 206 under ETag v1.</summary>
    private static HttpResponseMessage RestFrom(int from)
    {
        var response = new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = new ByteArrayContent(Body[from..]) };
        response.Content.Headers.ContentRange = new ContentRangeHeaderValue(from, Body.Length - 1, Body.Length);
        response.Headers.ETag = EntityTagHeaderValue.Parse("\"v1\"");
        return response;
    }

    [Test]
    public async Task Policy_HttpsLocation_PicksEachHttpsTry_GivenThePreviousTry()
    {
        // PXReprise's route (040): Globus HTTPS first; after a Globus try that fails without a damaged copy, one try
        // from ftp.pride; then Globus again.
        int globusRequests = 0;
        var handler = new StubHandler(request => request.RequestUri.Host == GlobusHost
            ? ++globusRequests == 1 ? new HttpResponseMessage(HttpStatusCode.ServiceUnavailable) : Download(Body)
            : new HttpResponseMessage(HttpStatusCode.BadGateway));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);
        var asked = new List<PrideDownloadTry>();
        var tries = new List<PrideDownloadTry>();

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir, Expected, new PrideDownloadPolicy
        {
            Tries = new SyncProgress<PrideDownloadTry>(tries.Add),
            HttpsLocation = (_, previous) =>
            {
                asked.Add(previous);
                return previous is { Host: GlobusHost, Outcome: PrideDownloadOutcome.Failed } ? null : GlobusUrl;
            }
        });

        string[] hosts = { GlobusHost, "ftp.pride.ebi.ac.uk", GlobusHost };
        Assert.That(handler.Requests.Select(r => r.RequestUri.Host), Is.EqualTo(hosts), "null is the file's own location");
        Assert.That(tries.Select(t => t.Host), Is.EqualTo(hosts));
        Assert.That(tries.Select(t => t.Outcome), Is.EqualTo(new[]
            { PrideDownloadOutcome.Failed, PrideDownloadOutcome.Failed, PrideDownloadOutcome.Accepted }));
        Assert.That(asked, Is.EqualTo(new[] { null, tries[0], tries[1] }));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
    }

    [Test]
    public async Task Policy_HttpsLocation_IsAskedOnlyForHttpsTries()
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientWith(Says(1, "Session Stop  (Error: Session data transfer timeout)"), handler);
        var asked = new List<PrideDownloadTry>();

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected, new PrideDownloadPolicy
        {
            Aspera = Options,
            HttpsLocation = (_, previous) => { asked.Add(previous); return null; }
        });

        Assert.That(_started, Has.Count.EqualTo(3));
        Assert.That(asked, Has.Count.EqualTo(1));
        Assert.That(asked[0].Route, Is.EqualTo(PrideDownloadRoute.Aspera));
        Assert.That(asked[0].Try, Is.EqualTo(3));
    }

    [TestCase(false)]
    [TestCase(true)]
    public async Task Policy_HttpsLocation_OnAnotherHost_ResumesAPartialCopyOnlyWhenItsMd5IsChecked(bool verifyMd5)
    {
        // The partial copy was started on Globus, under Globus's validator; the next try is on ftp.pride.
        int call = 0;
        var handler = new StubHandler(request => ++call == 1 ? DropsAfter(4)
            : request.Headers.Range != null ? RestFrom(4) : Download(Body));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir, Expected,
            new PrideDownloadPolicy { HttpsLocation = (_, previous) => previous == null ? GlobusUrl : null },
            verifyMd5: verifyMd5);

        HttpRequestMessage second = handler.Requests[1];
        Assert.That(second.RequestUri.Host, Is.EqualTo("ftp.pride.ebi.ac.uk"));
        Assert.That(second.Headers.IfRange, Is.Null, "another host's validator is never sent");
        Assert.That(second.Headers.Range?.ToString(), verifyMd5 ? Is.EqualTo("bytes=4-") : Is.Null);
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
    }

    [TestCase("http://globus.example.org/run1.raw")]
    [TestCase("ftp://globus.example.org/run1.raw")]
    [TestCase("globus.example.org/run1.raw")]
    public void Policy_HttpsLocation_ThatIsNotAnAbsoluteHttpsUrl_IsRefusedAtOnce(string location)
    {
        var handler = new StubHandler(_ => Download(Body));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);
        var tries = new List<PrideDownloadTry>();

        var e = Assert.ThrowsAsync<ArgumentException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir,
            Expected, new PrideDownloadPolicy { Tries = new SyncProgress<PrideDownloadTry>(tries.Add), HttpsLocation = (_, _) => location }));

        Assert.That(e.ParamName, Is.EqualTo("policy"));
        Assert.That(e.Message, Does.Not.Contain(GlobusHost), "the location may carry a credential");
        Assert.That(handler.Requests, Is.Empty);
        Assert.That(tries, Is.Empty);
    }

    [Test]
    public async Task Policy_Tries_RecordWhenTheyStarted_AndTheBytesTheyReceived()
    {
        // Aspera fails having written nothing, then writes a damaged 11-byte copy; HTTPS then brings the 10 bytes.
        int run = 0;
        using var client = ClientWith(psi => ++run == 1
                ? Says(1, "Session Stop  (Error: Session data transfer timeout)")(psi)
                : Writes(Encoding.ASCII.GetBytes("eleven byte"))(psi),
            new StubHandler(_ => Download(Body)));
        var tries = new List<PrideDownloadTry>();
        DateTime before = DateTime.UtcNow;

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected,
            new PrideDownloadPolicy { Aspera = Options, Tries = new SyncProgress<PrideDownloadTry>(tries.Add) });

        DateTime after = DateTime.UtcNow;
        Assert.That(tries.Select(t => (t.Route, t.Outcome, t.BytesReceived)), Is.EqualTo(new[]
        {
            (PrideDownloadRoute.Aspera, PrideDownloadOutcome.Failed, (long?)0),
            (PrideDownloadRoute.Aspera, PrideDownloadOutcome.Damaged, (long?)11),
            (PrideDownloadRoute.Https, PrideDownloadOutcome.Accepted, (long?)Body.Length)
        }));
        Assert.That(tries.Select(t => t.StartedUtc.Kind), Is.All.EqualTo(DateTimeKind.Utc));
        Assert.That(tries.Select(t => t.StartedUtc), Is.All.InRange(before, after));
        Assert.That(tries.Select(t => t.StartedUtc), Is.Ordered);
    }

    [Test]
    public async Task Policy_AcceptedAsperaCopy_ReceivedItsSize()
    {
        using var client = ClientWith(Writes(Body));
        var tries = new List<PrideDownloadTry>();

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir, Expected,
            new PrideDownloadPolicy { Aspera = Options, Tries = new SyncProgress<PrideDownloadTry>(tries.Add) });

        Assert.That(tries.Single().Outcome, Is.EqualTo(PrideDownloadOutcome.Accepted));
        Assert.That(tries.Single().BytesReceived, Is.EqualTo(Body.Length));
    }

    [Test]
    public async Task Policy_RelativeDestination_OverAspera_ReturnsTheFullPath()
    {
        string relative = RelativeTempDir();
        using var client = ClientWith(Writes(Body));

        string path = await client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), relative, Expected,
            new PrideDownloadPolicy { Aspera = Options, Attempts = 1 });

        Assert.That(path, Is.EqualTo(Destination), "the full path, as documented");
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
    }

    [Test]
    public void Policy_StalledAsperaTry_RecordsItsStart_AndTheBytesBeforeTheStall()
    {
        // ascp writes three bytes and then nothing, so the session is stopped once the stall window passes.
        var stall = TimeSpan.FromMilliseconds(300);
        using var client = ClientWith(psi =>
        {
            File.WriteAllBytes(Path.Combine(StagingOf(psi), "run1.raw"), Body[..3]);
            return new FakeAscp();
        });
        var reports = new List<(PrideDownloadTry Try, DateTime ReportedUtc)>();

        Assert.ThrowsAsync<PrideAsperaException>(() => client.DownloadFileWithPolicyAsync(MakeFile(Ftp, Aspera), _tempDir,
            Expected, new PrideDownloadPolicy
            {
                Aspera = new PrideAsperaOptions { AscpPath = Options.AscpPath, StallTimeout = stall },
                Attempts = 1,
                Tries = new SyncProgress<PrideDownloadTry>(t => reports.Add((t, DateTime.UtcNow)))
            }));

        (PrideDownloadTry stalled, DateTime reportedUtc) = reports.Single();
        Assert.That(stalled.AsperaFailure, Is.EqualTo(PrideAsperaFailure.Stalled));
        Assert.That(stalled.BytesReceived, Is.EqualTo(3));
        Assert.That(stalled.Elapsed, Is.GreaterThanOrEqualTo(stall));
        Assert.That(reportedUtc - stalled.StartedUtc, Is.GreaterThanOrEqualTo(stall), "when it started, not when it ended");
    }

    [Test]
    public async Task Policy_ResumedHttpsTry_CountsOnlyTheBytesItReceived()
    {
        int call = 0;
        var handler = new StubHandler(_ => ++call == 1 ? DropsAfter(4) : RestFrom(4));
        using var client = ClientWith(_ => throw new AssertionException("ascp must not run"), handler);
        var tries = new List<PrideDownloadTry>();

        await client.DownloadFileWithPolicyAsync(MakeFile(Ftp), _tempDir, Expected,
            new PrideDownloadPolicy { Tries = new SyncProgress<PrideDownloadTry>(tries.Add) });

        HttpRequestMessage second = handler.Requests[1];
        Assert.That(second.Headers.Range?.ToString(), Is.EqualTo("bytes=4-"), "the same host resumes");
        Assert.That(second.Headers.IfRange?.EntityTag?.Tag, Is.EqualTo("\"v1\""));
        Assert.That(tries.Select(t => t.Outcome), Is.EqualTo(new[] { PrideDownloadOutcome.Failed, PrideDownloadOutcome.Accepted }));
        Assert.That(tries.Select(t => t.BytesReceived), Is.EqualTo(new long?[] { 4, Body.Length - 4 }));
        Assert.That(File.ReadAllBytes(Destination), Is.EqualTo(Body));
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
