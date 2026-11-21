using System;
using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Runtime.ExceptionServices;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    public sealed partial class PrideArchiveClient
    {
        /// <summary>The manifest name of EBI's public Aspera key, embedded from <c>Resources/aspera_tokenauth_id_rsa</c>.</summary>
        internal const string AsperaKeyResourceName = "UsefulProteomicsDatabases.Resources.aspera_tokenauth_id_rsa";

        /// <summary>The folder, inside the destination directory, that holds each Aspera attempt's private folder.</summary>
        internal const string AsperaFolderName = ".aspera";

        /// <summary>ascp's port for EBI's server (EBI KB0011565).</summary>
        private const string AsperaPort = "33001";

        /// <summary>Starts ascp. Tests replace it to script ascp's behaviour without running it.</summary>
        internal Func<ProcessStartInfo, IAscpSession> StartAscp { get; init; } = AscpProcess.Start;

        /// <summary>The Aspera route's cap, start gap and pause, shared by the process. Tests give each client its own.</summary>
        internal AsperaRoute AsperaRouteState { get; init; } = AsperaRoute.Shared;

        /// <summary>
        /// Downloads <paramref name="file"/> as <paramref name="policy"/> says: over Aspera and HTTPS, retrying with
        /// growing waits, setting Aspera aside after a route failure, and reporting each try. Every copy is checked
        /// against PRIDE's checksum row before it is moved into place.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Opt-in: the other downloads make one Aspera attempt, or one HTTPS download with a fixed internal retry, and
        /// report what went wrong. This one decides what to do about it, the way a pipeline that downloads PRIDE for
        /// weeks needs to (the defaults, <see cref="PrideDownloadPolicy"/>, are PXReprise's).
        /// </para>
        /// <para>
        /// Each try goes over Aspera when <see cref="PrideDownloadPolicy.Aspera"/> is set, the file has an Aspera
        /// location, it has had fewer than <see cref="PrideDownloadPolicy.AsperaTries"/> Aspera tries, no Aspera copy of
        /// it came back damaged (when <see cref="PrideDownloadPolicy.HttpsAfterDamagedAsperaCopy"/>), and the route is
        /// not set aside; otherwise over HTTPS. A file with no HTTPS location goes over Aspera whatever those say, since
        /// it has no other route. An Aspera try waits for a session under the process's cap and start gap.
        /// </para>
        /// <para>
        /// An HTTPS try is one attempt, not the HTTPS download's own retries: this method's tries replace them, so the
        /// two budgets do not multiply. Between tries it keeps what the HTTPS download keeps (a <c>.partial</c> to
        /// resume after an outage) and never resumes an Aspera copy. An HTTPS try downloads from the file's own HTTPS
        /// location, or from where <see cref="PrideDownloadPolicy.HttpsLocation"/> sends it.
        /// </para>
        /// <para>
        /// Tried again: an outage on either route (<see cref="ServiceFailure.IsTransient"/>), any
        /// <see cref="PrideAsperaException"/>, and a copy that fails its check. Not tried again: a refusal such as 404,
        /// a broken contract, the caller's disk, a configuration error and cancellation, which are thrown at once.
        /// After the last try, its failure is thrown: an <see cref="HttpRequestException"/> (perhaps a
        /// <see cref="PrideAsperaException"/>) for an outage, an <see cref="MzLibException"/> for a copy that failed its
        /// check over HTTPS.
        /// </para>
        /// </remarks>
        /// <param name="file">The file to download.</param>
        /// <param name="destinationDirectory">The directory to write into; created if it does not exist. A relative path is taken from the current directory when the call starts.</param>
        /// <param name="expected">PRIDE's checksum row for this file. Its <see cref="PrideFileChecksum.FileName"/> must equal <paramref name="file"/>'s.</param>
        /// <param name="policy">The routes, tries, waits and pause.</param>
        /// <param name="overwrite">
        /// When true (the default) an existing destination file is replaced. When false, an existing destination file
        /// that passes the check is left untouched and nothing is downloaded.
        /// </param>
        /// <param name="verifyMd5">When true, the MD5 is checked as well as the size. False by default.</param>
        /// <param name="cancellationToken">Cancels the try in progress, or the wait between tries.</param>
        /// <returns>The full path of the written (or already-present and checked) file.</returns>
        /// <exception cref="ArgumentNullException">The file, <paramref name="expected"/> or <paramref name="policy"/> is null.</exception>
        /// <exception cref="ArgumentException">As for the other downloads; or a <paramref name="policy"/> setting is out of range; or <see cref="PrideDownloadPolicy.HttpsLocation"/> returned something other than an absolute <c>https</c> URL.</exception>
        /// <exception cref="NotSupportedException">The file has neither an HTTPS-reachable location nor, with <see cref="PrideDownloadPolicy.Aspera"/> set, an Aspera one.</exception>
        /// <exception cref="HttpRequestException">Every try failed, the last with an outage.</exception>
        /// <exception cref="MzLibException">Every try failed, the last over HTTPS with a copy that failed its check; or PRIDE broke its contract; or ascp could not be started.</exception>
        /// <exception cref="IOException">The caller's disk failed.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled.</exception>
        public async Task<string> DownloadFileWithPolicyAsync(PrideArchiveFile file, string destinationDirectory,
            PrideFileChecksum expected, PrideDownloadPolicy policy, bool overwrite = true, bool verifyMd5 = false,
            CancellationToken cancellationToken = default)
        {
            if (expected == null)
                throw new ArgumentNullException(nameof(expected));
            CheckPolicy(policy);
            string safeFileName = CheckDownloadTarget(file, destinationDirectory, expected);
            // Resolved once, so every try writes into the same folder however long the waits between them, and the
            // path returned is the full one.
            destinationDirectory = Path.GetFullPath(destinationDirectory);
            string destinationPath = Path.Combine(destinationDirectory, safeFileName);

            if (!overwrite && File.Exists(destinationPath)
                && await ChecksumMismatchAsync(destinationPath, expected, verifyMd5, cancellationToken).ConfigureAwait(false) == null)
                return destinationPath;

            bool hasAspera = policy.Aspera != null && file.TryGetAsperaSource(out _);
            bool hasHttps = file.TryGetHttpsDownloadUrl(out string httpsUrl);
            if (!hasAspera && !hasHttps)
                throw new NotSupportedException(policy.Aspera == null
                    ? $"PRIDE file '{safeFileName}' exposes no HTTPS-reachable download location (only Aspera, or none), and the policy sets no Aspera options."
                    : $"PRIDE file '{safeFileName}' has neither an HTTPS-reachable download location nor an Aspera location on {PrideArchiveExtensions.PrideAsperaHost}.");

            int asperaTries = 0;
            bool asperaDamaged = false;
            string lastHttpsHost = null;
            PrideDownloadTry previous = null;
            Exception last = null;
            for (int attempt = 1; attempt <= policy.Attempts; attempt++)
            {
                if (attempt > 1)
                    await RetryDelay(policy.Backoff * (attempt - 1), cancellationToken).ConfigureAwait(false);

                bool viaAspera = hasAspera && (!hasHttps
                    || (asperaTries < policy.AsperaTries
                        && !(asperaDamaged && policy.HttpsAfterDamagedAsperaCopy)
                        && !AsperaRouteState.IsPaused));
                PrideDownloadRoute route = viaAspera ? PrideDownloadRoute.Aspera : PrideDownloadRoute.Https;
                string url = viaAspera ? null : HttpsLocationFor(policy, file, previous) ?? httpsUrl;
                string host = viaAspera ? PrideArchiveExtensions.PrideAsperaHost : _engine.HostOf(url);
                if (!viaAspera)
                {
                    // A partial copy's validator is the host's that started it, and vouches for nothing on another
                    // host. Without it the copy is resumed only when its MD5 is checked, as for any partial copy.
                    if (lastHttpsHost != null && !string.Equals(host, lastHttpsHost, StringComparison.OrdinalIgnoreCase))
                        DownloadEngine.DeleteQuietly(destinationPath + ".partial" + ValidatorSuffix);
                    lastHttpsHost = host;
                }

                DateTime? startedUtc = null;
                Stopwatch clock = null;
                var bytes = new TryBytes();
                try
                {
                    if (viaAspera)
                    {
                        asperaTries++;
                        using (await AsperaRouteState.EnterAsync(policy.MaxConcurrentAsperaSessions, policy.AsperaStartGap,
                            RetryDelay, cancellationToken).ConfigureAwait(false))
                        {
                            startedUtc = DateTime.UtcNow;
                            clock = Stopwatch.StartNew();
                            await DownloadFileViaAsperaCoreAsync(file, destinationDirectory, expected, policy.Aspera,
                                overwrite: true, verifyMd5, cancellationToken, bytes).ConfigureAwait(false);
                        }
                    }
                    else
                    {
                        startedUtc = DateTime.UtcNow;
                        clock = Stopwatch.StartNew();
                        await DownloadFileCoreAsync(file, destinationDirectory, expected, overwrite: true, verifyMd5,
                            cancellationToken, singleAttempt: true, httpsUrl: url, bytes: bytes).ConfigureAwait(false);
                    }

                    ReportTry(policy, safeFileName, attempt, route, host, expected, startedUtc, clock, bytes,
                        PrideDownloadOutcome.Accepted, null);
                    return destinationPath;
                }
                catch (Exception e) when (IsWorthAnotherTry(e, out PrideDownloadOutcome outcome))
                {
                    previous = ReportTry(policy, safeFileName, attempt, route, host, expected, startedUtc, clock, bytes,
                        outcome, e);
                    if (e is PrideAsperaException { Failure: PrideAsperaFailure.LoginRefused or PrideAsperaFailure.CannotConnect })
                        AsperaRouteState.Pause(policy.AsperaPauseAfterRouteFailure);
                    if (e is PrideAsperaException { Failure: PrideAsperaFailure.DamagedCopy })
                        asperaDamaged = true;
                    last = e;
                }
                catch (Exception e) when (e is not OperationCanceledException)
                {
                    ReportTry(policy, safeFileName, attempt, route, host, expected, startedUtc, clock, bytes,
                        PrideDownloadOutcome.Failed, e);
                    throw;
                }
            }

            ExceptionDispatchInfo.Throw(last);
            return null; // unreachable
        }

        /// <summary>The bytes one try received, filled in as it ends; null when they could not be counted.</summary>
        private sealed class TryBytes
        {
            public long? Received { get; set; }
        }

        /// <summary>
        /// Where <see cref="PrideDownloadPolicy.HttpsLocation"/> sends the next HTTPS try, or null for the file's own
        /// location. Neither the URL nor anything else from it goes into the message: it may carry a credential.
        /// </summary>
        private static string HttpsLocationFor(PrideDownloadPolicy policy, PrideArchiveFile file, PrideDownloadTry previous)
        {
            string url = policy.HttpsLocation?.Invoke(file, previous);
            if (url != null && !(Uri.TryCreate(url, UriKind.Absolute, out Uri uri) && uri.Scheme == Uri.UriSchemeHttps))
                throw new ArgumentException(
                    "PrideDownloadPolicy.HttpsLocation returned a location that is not an absolute https URL.", nameof(policy));
            return url;
        }

        /// <summary>
        /// Whether a try's failure is one another try may mend, and how the try ended: an outage on either route, any
        /// Aspera failure, or a copy that failed its check over HTTPS (the HTTPS download reports that, after its one
        /// attempt, as an <see cref="MzLibException"/> around a <see cref="ContentMismatchException"/>).
        /// </summary>
        private static bool IsWorthAnotherTry(Exception e, out PrideDownloadOutcome outcome)
        {
            outcome = PrideDownloadOutcome.Failed;
            switch (e)
            {
                case PrideAsperaException asperaFailure:
                    if (asperaFailure.Failure == PrideAsperaFailure.DamagedCopy)
                        outcome = PrideDownloadOutcome.Damaged;
                    return true;
                case MzLibException { InnerException: ContentMismatchException }:
                    outcome = PrideDownloadOutcome.Damaged;
                    return true;
                case HttpRequestException:
                    return ServiceFailure.IsTransient(e);
                default:
                    return false;
            }
        }

        private static PrideDownloadTry ReportTry(PrideDownloadPolicy policy, string fileName, int attempt,
            PrideDownloadRoute route, string host, PrideFileChecksum expected, DateTime? startedUtc, Stopwatch clock,
            TryBytes bytes, PrideDownloadOutcome outcome, Exception error)
        {
            var record = new PrideDownloadTry(fileName, attempt, route, host, expected.SizeBytes,
                startedUtc ?? DateTime.UtcNow, clock?.Elapsed ?? TimeSpan.Zero, bytes.Received, outcome, error);
            policy.Tries?.Report(record);
            return record;
        }

        private static void CheckPolicy(PrideDownloadPolicy policy)
        {
            if (policy == null)
                throw new ArgumentNullException(nameof(policy));
            if (policy.Aspera != null)
                CheckAsperaOptions(policy.Aspera);
            if (policy.Attempts < 1)
                throw new ArgumentException("PrideDownloadPolicy.Attempts must be at least 1.", nameof(policy));
            if (policy.AsperaTries < 0)
                throw new ArgumentException("PrideDownloadPolicy.AsperaTries must not be negative.", nameof(policy));
            if (policy.MaxConcurrentAsperaSessions < 1)
                throw new ArgumentException("PrideDownloadPolicy.MaxConcurrentAsperaSessions must be at least 1.", nameof(policy));
            if (policy.Backoff < TimeSpan.Zero || policy.AsperaStartGap < TimeSpan.Zero || policy.AsperaPauseAfterRouteFailure < TimeSpan.Zero)
                throw new ArgumentException(
                    "PrideDownloadPolicy.Backoff, AsperaStartGap and AsperaPauseAfterRouteFailure must not be negative.", nameof(policy));
        }

        /// <summary>
        /// Downloads <paramref name="file"/> from EBI's Aspera server with IBM Aspera's <c>ascp</c>, once, and checks it
        /// against PRIDE's checksum row before moving it into <paramref name="destinationDirectory"/>.
        /// </summary>
        /// <remarks>
        /// <para>
        /// One attempt, and nothing else: no retry, no fall-back to HTTPS, no limit on how many run at once. A failure
        /// is reported with its kind (<see cref="PrideAsperaException.Failure"/>) so the caller can decide what to do
        /// next. <see cref="DownloadFileWithPolicyAsync"/> is the caller that does those things, when asked.
        /// </para>
        /// <para>
        /// The command is <c>ascp -Q -l &lt;MaxRate&gt; -P 33001 -i &lt;key&gt; &lt;source&gt; &lt;folder&gt;</c>, EBI's (KB0011565):
        /// encryption stays on, because <c>-T</c> gave a damaged file. ascp writes into a private folder,
        /// <c>.aspera/&lt;guid&gt;/</c> inside <paramref name="destinationDirectory"/>, which is always deleted afterwards. An
        /// attempt never resumes an earlier one: resuming after a killed session gave damaged files (PXReprise,
        /// 2026-10-08).
        /// </para>
        /// <para>
        /// The copy is checked as <see cref="DownloadFileAsync(PrideArchiveFile, string, PrideFileChecksum, bool, bool, CancellationToken)"/>
        /// checks an HTTPS download: its size, then its MD5 when <paramref name="verifyMd5"/> is set, else a scan for
        /// EBI's error text, which EBI has written into files served over Aspera as well as HTTPS. With the MD5, a copy
        /// holding EBI's exact opening (<c>&lt;Error&gt;&lt;Code&gt;ConnectionClosedException&lt;/Code&gt;</c>) is refused even
        /// when its MD5 matches, as an HTTPS download stopped at that opening is. ascp writes the file itself, so nothing
        /// can be checked while it arrives: the whole file is read once at the end. Only a copy that
        /// passes is moved to its destination; it then replaces any <c>.partial</c> an HTTPS download left for the same
        /// file.
        /// </para>
        /// <para>
        /// A session whose file receives no bytes for <see cref="PrideAsperaOptions.StallTimeout"/> is stopped, with every
        /// process ascp started. Cancelling stops it the same way.
        /// </para>
        /// </remarks>
        /// <param name="file">The file to download. It must have an Aspera location on <see cref="PrideArchiveExtensions.PrideAsperaHost"/>.</param>
        /// <param name="destinationDirectory">The directory to write into; created if it does not exist. A relative path is taken from the current directory when the call starts.</param>
        /// <param name="expected">PRIDE's checksum row for this file. Its <see cref="PrideFileChecksum.FileName"/> must equal <paramref name="file"/>'s.</param>
        /// <param name="aspera">Where ascp is, and how to run it.</param>
        /// <param name="overwrite">
        /// When true (the default) an existing destination file is replaced. When false, an existing destination file
        /// that passes the check is left untouched and ascp is not run.
        /// </param>
        /// <param name="verifyMd5">When true, the MD5 is checked as well as the size. False by default.</param>
        /// <param name="cancellationToken">Stops ascp, or the check.</param>
        /// <returns>The full path of the written (or already-present and checked) file.</returns>
        /// <exception cref="ArgumentNullException">The file, <paramref name="expected"/> or <paramref name="aspera"/> is null.</exception>
        /// <exception cref="ArgumentException">As for the HTTPS download, or <paramref name="aspera"/> has no <see cref="PrideAsperaOptions.AscpPath"/> or an unreadable <see cref="PrideAsperaOptions.MaxRate"/>, or a blank <see cref="PrideAsperaOptions.KeyPath"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><see cref="PrideAsperaOptions.StallTimeout"/> is not positive.</exception>
        /// <exception cref="NotSupportedException">The file has no Aspera location on <see cref="PrideArchiveExtensions.PrideAsperaHost"/>.</exception>
        /// <exception cref="PrideAsperaException">The attempt failed; <see cref="PrideAsperaException.Failure"/> says how.</exception>
        /// <exception cref="MzLibException">ascp could not be started, or <see cref="PrideAsperaOptions.KeyPath"/> does not exist: a configuration error, not an outage.</exception>
        /// <exception cref="IOException">ascp ran out of disk space, or the copy could not be moved into place.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled.</exception>
        public Task<string> DownloadFileViaAsperaAsync(PrideArchiveFile file, string destinationDirectory,
            PrideFileChecksum expected, PrideAsperaOptions aspera, bool overwrite = true, bool verifyMd5 = false,
            CancellationToken cancellationToken = default) =>
            DownloadFileViaAsperaCoreAsync(file, destinationDirectory, expected, aspera, overwrite, verifyMd5,
                cancellationToken, bytes: null);

        /// <param name="bytes">Receives the bytes in ascp's folder when it ended; null when nobody counts them.</param>
        private async Task<string> DownloadFileViaAsperaCoreAsync(PrideArchiveFile file, string destinationDirectory,
            PrideFileChecksum expected, PrideAsperaOptions aspera, bool overwrite, bool verifyMd5,
            CancellationToken cancellationToken, TryBytes bytes)
        {
            if (expected == null)
                throw new ArgumentNullException(nameof(expected));
            CheckAsperaOptions(aspera);
            string safeFileName = CheckDownloadTarget(file, destinationDirectory, expected);
            // ascp runs in its staging folder (WorkingDirectory below), where a relative path names another file. So the
            // folder, and the key below, are resolved here, against this process's directory, before ascp is given them.
            destinationDirectory = Path.GetFullPath(destinationDirectory);
            string destinationPath = Path.Combine(destinationDirectory, safeFileName);

            if (!overwrite && File.Exists(destinationPath)
                && await ChecksumMismatchAsync(destinationPath, expected, verifyMd5, cancellationToken).ConfigureAwait(false) == null)
                return destinationPath;

            if (!file.TryGetAsperaSource(out string source))
                throw new NotSupportedException(
                    $"PRIDE file '{safeFileName}' has no Aspera location on {PrideArchiveExtensions.PrideAsperaHost}.");

            // ascp names the copy after the source; TryGetAsperaSource guarantees that is a plain leaf name.
            string sourceName = source.Substring(source.LastIndexOf('/') + 1);
            string described = $"'{safeFileName}' from {PrideArchiveExtensions.PrideAsperaHost}";

            Directory.CreateDirectory(destinationDirectory);
            string asperaRoot = Path.Combine(destinationDirectory, AsperaFolderName);
            string attemptDirectory = Path.Combine(asperaRoot, Guid.NewGuid().ToString("N"));
            string stagingDirectory = Path.Combine(attemptDirectory, "files");
            Directory.CreateDirectory(stagingDirectory);

            try
            {
                string keyPath = aspera.KeyPath != null ? Path.GetFullPath(aspera.KeyPath) : WriteBundledAsperaKey(attemptDirectory);
                if (!File.Exists(keyPath))
                    throw new MzLibException($"The Aspera key file '{keyPath}' does not exist.");

                var startInfo = new ProcessStartInfo(aspera.AscpPath)
                {
                    UseShellExecute = false,
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                    CreateNoWindow = true,
                    WorkingDirectory = stagingDirectory
                };
                foreach (string argument in new[] { "-Q", "-l", aspera.MaxRate, "-P", AsperaPort, "-i", keyPath, source, stagingDirectory })
                    startInfo.ArgumentList.Add(argument);

                (int exitCode, string said) = await RunAscpAsync(startInfo, stagingDirectory, described, aspera, cancellationToken).ConfigureAwait(false);
                if (bytes != null)
                    bytes.Received = Counted(StagedBytes(stagingDirectory));
                said = Scrub(said, source, attemptDirectory, keyPath);

                if (exitCode != 0)
                    throw AscpExitFailure(described, exitCode, said);

                string stagedPath = Path.Combine(stagingDirectory, sourceName);
                if (!File.Exists(stagedPath))
                    throw new PrideAsperaException(
                        $"The PRIDE Aspera download of {described} failed: ascp exited with 0 but wrote no file.",
                        PrideAsperaFailure.AscpFailed, exitCode,
                        new IOException($"ascp exited with 0 and wrote no '{sourceName}'. Its last line: {said}"));

                // ascp wrote the file, so it was never seen as it arrived: read it once, whole.
                using (var digest = new DownloadDigest(verifyMd5))
                {
                    digest.Invalidate();
                    string mismatch = await ContentMismatchAsync(stagedPath, expected, verifyMd5, digest, cancellationToken).ConfigureAwait(false);

                    // Over HTTPS, EBI's exact opening stops the transfer as it arrives, so a file holding it is refused
                    // even when its MD5 would match (W1b). ascp cannot be stopped there; the re-read keeps the same rule.
                    if (mismatch == null && verifyMd5 && digest.ServerErrorAt >= 0)
                        mismatch = ServerErrorMismatch(digest.ServerErrorAt);
                    if (mismatch != null)
                        throw new PrideAsperaException($"The PRIDE Aspera download of {described} {mismatch}.",
                            PrideAsperaFailure.DamagedCopy, exitCode, new InvalidDataException($"The copy {mismatch}."));
                }

                File.Move(stagedPath, destinationPath, overwrite: true);

                // An HTTPS download of the same file that stopped part way is now out of date.
                string partialPath = destinationPath + ".partial";
                DownloadEngine.DeleteQuietly(partialPath);
                DownloadEngine.DeleteQuietly(partialPath + ValidatorSuffix);
                return destinationPath;
            }
            finally
            {
                // ascp stalled, could not start, or was cancelled: what it wrote is still in its folder, until now.
                if (bytes != null && bytes.Received == null)
                    bytes.Received = Counted(StagedBytes(stagingDirectory));
                DeleteDirectoryQuietly(attemptDirectory, recursive: true);
                DeleteDirectoryQuietly(asperaRoot, recursive: false);
            }
        }

        /// <summary><see cref="StagedBytes"/>'s count, or null where it could not read the folder.</summary>
        private static long? Counted(long stagedBytes) => stagedBytes < 0 ? null : stagedBytes;

        /// <summary>The checks on <see cref="PrideAsperaOptions"/> that need no file system.</summary>
        private static void CheckAsperaOptions(PrideAsperaOptions aspera)
        {
            if (aspera == null)
                throw new ArgumentNullException(nameof(aspera));
            if (string.IsNullOrWhiteSpace(aspera.AscpPath))
                throw new ArgumentException("The path of ascp is required (PrideAsperaOptions.AscpPath).", nameof(aspera));
            if (aspera.KeyPath != null && string.IsNullOrWhiteSpace(aspera.KeyPath))
                throw new ArgumentException(
                    "PrideAsperaOptions.KeyPath is blank: leave it null for the bundled key, or give the key file's path.",
                    nameof(aspera));
            if (!IsAscpRate(aspera.MaxRate))
                throw new ArgumentException(
                    $"'{aspera.MaxRate}' is not a rate ascp takes: digits, optionally followed by k, m or g (PrideAsperaOptions.MaxRate).",
                    nameof(aspera));
            if (aspera.StallTimeout <= TimeSpan.Zero)
                throw new ArgumentOutOfRangeException(nameof(aspera), aspera.StallTimeout,
                    "PrideAsperaOptions.StallTimeout must be positive.");
        }

        private static bool IsAscpRate(string rate)
        {
            if (string.IsNullOrEmpty(rate))
                return false;
            string digits = "kKmMgG".Contains(rate[^1]) ? rate[..^1] : rate;
            return digits.Length > 0 && digits.All(char.IsAsciiDigit);
        }

        /// <summary>
        /// Writes EBI's public Aspera key, which mzLib embeds, into <paramref name="attemptDirectory"/> for ascp to read,
        /// readable by its owner only where the OS has such modes. Deleted with the folder.
        /// </summary>
        private static string WriteBundledAsperaKey(string attemptDirectory)
        {
            string keyPath = Path.Combine(attemptDirectory, "aspera_tokenauth_id_rsa");
            using (Stream key = typeof(PrideArchiveClient).Assembly.GetManifestResourceStream(AsperaKeyResourceName)
                ?? throw new MzLibException(
                    $"The embedded resource '{AsperaKeyResourceName}' is missing; check the EmbeddedResource entry in UsefulProteomicsDatabases.csproj."))
            using (var file = new FileStream(keyPath, FileMode.CreateNew, FileAccess.Write))
            {
                key.CopyTo(file);
            }

            if (!OperatingSystem.IsWindows())
                File.SetUnixFileMode(keyPath, UnixFileMode.UserRead | UnixFileMode.UserWrite);
            return keyPath;
        }

        /// <summary>
        /// Runs ascp until it exits, watching its folder: no new bytes for <see cref="PrideAsperaOptions.StallTimeout"/>
        /// stops it as <see cref="PrideAsperaFailure.Stalled"/>, and cancellation stops it too.
        /// </summary>
        private async Task<(int ExitCode, string LastLine)> RunAscpAsync(ProcessStartInfo startInfo, string stagingDirectory,
            string described, PrideAsperaOptions aspera, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();

            IAscpSession session;
            try
            {
                session = StartAscp(startInfo);
            }
            catch (Exception e) when (e is Win32Exception or InvalidOperationException or IOException)
            {
                throw new MzLibException($"ascp could not be started from '{aspera.AscpPath}': {e.Message}", e);
            }

            using (session)
            {
                // Polled often enough to notice a stall within a quarter of its window; every 15 s at most.
                TimeSpan poll = TimeSpan.FromTicks(Math.Min(TimeSpan.FromSeconds(15).Ticks, Math.Max(aspera.StallTimeout.Ticks / 4, TimeSpan.FromMilliseconds(10).Ticks)));
                long seen = StagedBytes(stagingDirectory);
                var quiet = Stopwatch.StartNew();
                while (true)
                {
                    Task tick = Task.Delay(poll, cancellationToken);
                    if (await Task.WhenAny(session.Exited, tick).ConfigureAwait(false) == session.Exited)
                        return (await session.Exited.ConfigureAwait(false), session.LastLine);

                    if (cancellationToken.IsCancellationRequested)
                    {
                        await StopAsync(session).ConfigureAwait(false);
                        cancellationToken.ThrowIfCancellationRequested();
                    }

                    long now = StagedBytes(stagingDirectory);
                    if (now != seen)
                    {
                        seen = now;
                        quiet.Restart();
                        continue;
                    }

                    if (quiet.Elapsed >= aspera.StallTimeout)
                    {
                        await StopAsync(session).ConfigureAwait(false);
                        string stall = $"no bytes for {aspera.StallTimeout} ({Math.Max(now, 0)} bytes in all)";
                        throw new PrideAsperaException(
                            $"The PRIDE Aspera download of {described} received {stall}, so ascp was stopped.",
                            PrideAsperaFailure.Stalled, null, new TimeoutException($"ascp received {stall}."));
                    }
                }
            }
        }

        /// <summary>Stops ascp and waits up to 30 s for it to go, so its files can be deleted.</summary>
        private static async Task StopAsync(IAscpSession session)
        {
            session.Kill();
            try
            {
                await session.Exited.WaitAsync(TimeSpan.FromSeconds(30)).ConfigureAwait(false);
            }
            catch (TimeoutException)
            {
                // The folder is deleted best-effort; a process that will not die is not worth blocking the caller on.
            }
        }

        /// <summary>The bytes in ascp's folder so far: it names its in-progress file itself. -1 when it cannot be read.</summary>
        private static long StagedBytes(string stagingDirectory)
        {
            try
            {
                return Directory.EnumerateFiles(stagingDirectory, "*", SearchOption.AllDirectories).Sum(f => new FileInfo(f).Length);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                return -1;
            }
        }

        /// <summary>
        /// The exception for ascp exiting with <paramref name="exitCode"/>, by what it said last (the phrases are ascp's
        /// own, as PXReprise's batch logged them from 2026-10-08 to 2026-11-15).
        /// </summary>
        internal static Exception AscpExitFailure(string described, int exitCode, string said)
        {
            string ascpSaid = $"ascp exited with {exitCode}: {said}";
            if (said.Contains("No space left on disk", StringComparison.OrdinalIgnoreCase))
                return new IOException($"The PRIDE Aspera download of {described} ran out of disk space ({ascpSaid}).");

            PrideAsperaFailure failure = ClassifyAscpFailure(said);
            string what = failure switch
            {
                PrideAsperaFailure.LoginRefused => "was refused: EBI did not accept the key",
                PrideAsperaFailure.CannotConnect => "could not connect to the server",
                _ => "failed"
            };
            return new PrideAsperaException($"The PRIDE Aspera download of {described} {what} ({ascpSaid}).",
                failure, exitCode, new IOException(ascpSaid));
        }

        /// <summary>What ascp's last line says went wrong.</summary>
        internal static PrideAsperaFailure ClassifyAscpFailure(string said)
        {
            if (said.Contains("failed to authenticate", StringComparison.OrdinalIgnoreCase))
                return PrideAsperaFailure.LoginRefused;
            if (said.Contains("Unable to connect via SSH", StringComparison.OrdinalIgnoreCase)
                || said.Contains("Failed to open TCP connection", StringComparison.OrdinalIgnoreCase)
                || said.Contains("unable to connect to server", StringComparison.OrdinalIgnoreCase)
                || said.Contains("Session initiation failed", StringComparison.OrdinalIgnoreCase))
                return PrideAsperaFailure.CannotConnect;
            return PrideAsperaFailure.AscpFailed;
        }

        /// <summary>
        /// ascp's line with the source and local paths replaced, so a message names neither (the rule
        /// <see cref="DownloadEngine.DescribeDownload"/> keeps for URLs).
        /// </summary>
        private static string Scrub(string said, string source, string attemptDirectory, string keyPath)
        {
            if (string.IsNullOrEmpty(said))
                return "(no output)";
            string remotePath = source.Substring(source.IndexOf(':') + 1);
            return said
                .Replace(source, "<source>", StringComparison.OrdinalIgnoreCase)
                .Replace(remotePath, "<source>", StringComparison.OrdinalIgnoreCase)
                .Replace(remotePath.TrimStart('/'), "<source>", StringComparison.OrdinalIgnoreCase)
                .Replace(keyPath, "<key>", StringComparison.OrdinalIgnoreCase)
                .Replace(attemptDirectory, "<folder>", StringComparison.OrdinalIgnoreCase);
        }

        private static void DeleteDirectoryQuietly(string path, bool recursive)
        {
            try
            {
                if (Directory.Exists(path) && (recursive || !Directory.EnumerateFileSystemEntries(path).Any()))
                    Directory.Delete(path, recursive);
            }
            catch (Exception e) when (e is IOException or UnauthorizedAccessException)
            {
                // Best effort: a leftover folder costs disk, not correctness, and the next attempt uses a new one.
            }
        }
    }
}
