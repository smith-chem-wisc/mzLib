using System;
using System.IO;
using System.Net;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Threading;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The transport a repository client downloads through: it retries transient failures, reports
    /// <see cref="HttpClient"/>'s own timeout and a stalled or dropped body as the transport failures they are, and
    /// fetches a file into a <c>.partial</c> that a later attempt can resume. Moved out of
    /// <see cref="PrideArchiveClient"/> unchanged, so that the MassIVE client can share it instead of copying it.
    /// </summary>
    /// <remarks>
    /// Repository-neutral by construction. The client supplies its service's name, which every message names
    /// ("PRIDE did not respond ..."), and any status its own hosts use for a transient failure on top of the
    /// shared verdict (PRIDE: 403 on its FTP host). What a download's bytes are checked against stays with the
    /// client too, behind <see cref="IDownloadInspector"/>.
    /// <para>
    /// The settings are properties read when a request is made, never copied at construction: a client's object
    /// initializer runs AFTER its constructor, so an engine that copied them there would keep the defaults.
    /// </para>
    /// </remarks>
    internal sealed class DownloadEngine
    {
        /// <summary>
        /// The wait before each retry of a transient failure: 5 s, 20 s, then 60 s, so a failure that persists
        /// gives up after about 85 s of backoff. See <see cref="WithRetryAsync{T}"/> for what counts as transient.
        /// </summary>
        /// <remarks>
        /// Deliberately internal: these waits stay fixed. A caller that wants its own tries, waits and routes opts in
        /// to <see cref="PrideDownloadPolicy"/>, whose tries replace these instead of multiplying them. A knob here can
        /// be published later, but a published one cannot be taken back.
        /// </remarks>
        internal static readonly TimeSpan[] RetryBackoff =
            { TimeSpan.FromSeconds(5), TimeSpan.FromSeconds(20), TimeSpan.FromSeconds(60) };

        /// <summary>The longest <c>Retry-After</c> a 429 or 503 is allowed to make a retry wait.</summary>
        private static readonly TimeSpan MaxRetryAfter = TimeSpan.FromSeconds(60);

        private const string RetryAfterKey = "DownloadEngine.RetryAfter";

        private readonly HttpClient _httpClient;
        private readonly string _serviceName;
        private readonly Func<HttpStatusCode, string, bool> _alsoTransient;

        /// <param name="httpClient">The client's <see cref="HttpClient"/>; the engine never disposes it.</param>
        /// <param name="serviceName">The service as messages name it, e.g. "PRIDE".</param>
        /// <param name="alsoTransient">
        /// Whether a status the shared verdict calls permanent is transient on a given host; null when the service
        /// has no such rule. See <see cref="IsTransient"/>.
        /// </param>
        internal DownloadEngine(HttpClient httpClient, string serviceName, Func<HttpStatusCode, string, bool> alsoTransient = null)
        {
            _httpClient = httpClient ?? throw new ArgumentNullException(nameof(httpClient));
            _serviceName = serviceName ?? throw new ArgumentNullException(nameof(serviceName));
            _alsoTransient = alsoTransient;
        }

        /// <summary>How many times a transient failure is retried after the first attempt. Tests set 0 to pin a single attempt.</summary>
        internal int MaxRetries { get; set; } = RetryBackoff.Length;

        /// <summary>Waits out a backoff. Tests replace it to record the schedule instead of sleeping through it.</summary>
        internal Func<TimeSpan, CancellationToken, Task> RetryDelay { get; set; } = Task.Delay;

        /// <summary>
        /// How long a download waits for the NEXT bytes of a response body before abandoning the transfer. Each
        /// read gets a fresh window, so this bounds silence, not duration. Why two minutes, and why an inactivity
        /// deadline: <see cref="PrideArchiveClient.BodyStallTimeout"/>.
        /// </summary>
        internal TimeSpan BodyStallTimeout { get; set; } = TimeSpan.FromMinutes(2);

        /// <summary>
        /// Sends a GET, reporting <see cref="HttpClient.Timeout"/> expiring as the transport failure it is.
        /// </summary>
        /// <remarks>
        /// When the client's own timeout fires, <see cref="HttpClient"/> throws a
        /// <see cref="TaskCanceledException"/> -- the same type a caller's cancellation produces -- so a service
        /// that never answers escaped the documented contract (and <c>ExternalServiceTestHelper</c>, which
        /// skips only on transport failures). HttpClient marks its own timeout with an inner
        /// <see cref="TimeoutException"/>, so only that is converted; a caller's cancellation, through the token or
        /// through <see cref="HttpClient.CancelPendingRequests"/> on an injected client, is not.
        /// </remarks>
        /// <param name="requestUri">The URI to fetch.</param>
        /// <param name="completionOption">When the returned task completes.</param>
        /// <param name="cancellationToken">The caller's token; its cancellation is never converted.</param>
        /// <param name="described">How the request is named in the message; defaults to the URI itself.</param>
        internal Task<HttpResponseMessage> GetAsync(string requestUri, HttpCompletionOption completionOption,
            CancellationToken cancellationToken, string described = null) =>
            SendAsync(new HttpRequestMessage(HttpMethod.Get, requestUri), completionOption, cancellationToken,
                described ?? $"'{requestUri}'");

        /// <summary>
        /// <see cref="GetAsync"/> for a request that carries headers of its own (a resumed download's
        /// <c>Range</c>). Takes ownership of <paramref name="request"/>.
        /// </summary>
        internal async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, HttpCompletionOption completionOption,
            CancellationToken cancellationToken, string described)
        {
            using (request)
            {
                try
                {
                    return await _httpClient.SendAsync(request, completionOption, cancellationToken).ConfigureAwait(false);
                }
                catch (TaskCanceledException e) when (e.InnerException is TimeoutException)
                {
                    throw new HttpRequestException(
                        $"{_serviceName} did not respond within {_httpClient.Timeout} for {described}.", e);
                }
            }
        }

        /// <summary>
        /// Fetches <paramref name="requestUri"/> in full, retrying transient failures, and returns the response
        /// only once it has a success status. Any other status throws, worded by <paramref name="describeFailure"/>.
        /// </summary>
        internal Task<HttpResponseMessage> GetSuccessAsync(string requestUri,
            Func<HttpResponseMessage, string> describeFailure, CancellationToken cancellationToken) =>
            WithRetryAsync(async () =>
            {
                HttpResponseMessage response = await GetAsync(requestUri, HttpCompletionOption.ResponseContentRead, cancellationToken).ConfigureAwait(false);
                if (response.IsSuccessStatusCode)
                    return response;

                using (response)
                    throw StatusFailure(response, describeFailure(response));
            }, HostOf(requestUri), cancellationToken);

        /// <summary>
        /// The exception for a non-success status: an <see cref="HttpRequestException"/> carrying the status,
        /// plus the server's <c>Retry-After</c> when it sent one, for <see cref="WithRetryAsync{T}"/> to honour.
        /// </summary>
        internal static HttpRequestException StatusFailure(HttpResponseMessage response, string message)
        {
            var exception = new HttpRequestException(message, null, response.StatusCode);
            TimeSpan? retryAfter = response.Headers.RetryAfter?.Delta
                ?? (response.Headers.RetryAfter?.Date is DateTimeOffset date ? date - DateTimeOffset.UtcNow : null);
            if (retryAfter.HasValue)
                exception.Data[RetryAfterKey] = retryAfter.Value;
            return exception;
        }

        /// <summary>
        /// Runs <paramref name="attempt"/>, retrying it up to <see cref="MaxRetries"/> times while it fails
        /// transiently, and lets the last failure propagate unchanged: same type, same
        /// <see cref="HttpRequestException.StatusCode"/>, same message. A transient failure it gives up on is
        /// marked, so <see cref="ServiceFailure.IsTransient"/> can tell the caller that it is worth asking again later.
        /// </summary>
        /// <remarks>
        /// PR-A made every transport failure an <see cref="HttpRequestException"/>, so the decision is made on
        /// the exception alone (and the host, see <see cref="IsTransient"/>).
        /// <para>
        /// Never retried: an <c>MzLibException</c> (the service answered, and broke its contract, which a
        /// second ask will not mend), a write-side <see cref="IOException"/> (the caller's disk) and the
        /// caller's cancellation, which also ends a backoff in progress.
        /// </para>
        /// <para>
        /// The mark is made here, at the one place the verdict is reached, because it cannot be worked out again
        /// from the exception afterwards: a 403 is transient on PRIDE's FTP host only, and the exception does not
        /// carry its host.
        /// </para>
        /// </remarks>
        /// <param name="attempt">One complete attempt, from sending the request to checking the status.</param>
        /// <param name="host">The host the request goes to; decides whether a 403 is transient.</param>
        /// <param name="cancellationToken">Cancels the attempt and any backoff.</param>
        /// <param name="maxRetries">
        /// This call's retry limit in place of <see cref="MaxRetries"/>; null (the default) uses it. A caller that runs
        /// its own retries passes 0, so the two budgets do not multiply.
        /// </param>
        internal async Task<T> WithRetryAsync<T>(Func<Task<T>> attempt, string host, CancellationToken cancellationToken,
            int? maxRetries = null)
        {
            int limit = maxRetries ?? MaxRetries;
            for (int retry = 0; ; retry++)
            {
                try
                {
                    return await attempt().ConfigureAwait(false);
                }
                catch (HttpRequestException e) when (retry < limit && IsTransient(e, host))
                {
                    await RetryDelay(BackoffFor(retry, e), cancellationToken).ConfigureAwait(false);
                }
                catch (HttpRequestException e) when (IsTransient(e, host))
                {
                    ServiceFailure.MarkTransient(e);
                    throw;
                }
            }
        }

        /// <summary>Whether <see cref="WithRetryAsync{T}"/> would retry <paramref name="exception"/> against <paramref name="host"/>.</summary>
        /// <remarks>
        /// One with no status is a transport failure (a timeout, a dropped or stalled body) and is always transient.
        /// One with a status is transient for 408, 429 and every 5xx, and for whatever the service's own rule adds
        /// (PRIDE: 403 on its FTP host only, because EBI answers rate limiting there with 403 rather than 429).
        /// </remarks>
        internal bool IsTransient(HttpRequestException exception, string host)
        {
            if (exception.StatusCode is not HttpStatusCode status)
                return true;

            int code = (int)status;
            return status == HttpStatusCode.RequestTimeout
                || status == HttpStatusCode.TooManyRequests
                || code >= 500
                || (_alsoTransient?.Invoke(status, host) ?? false);
        }

        /// <summary>
        /// The wait before retry <paramref name="retry"/> (zero-based): the server's <c>Retry-After</c> on a
        /// 429 or 503, capped at <see cref="MaxRetryAfter"/>, else the fixed backoff.
        /// </summary>
        private static TimeSpan BackoffFor(int retry, HttpRequestException exception)
        {
            if (exception.StatusCode is HttpStatusCode.TooManyRequests or HttpStatusCode.ServiceUnavailable
                && exception.Data[RetryAfterKey] is TimeSpan retryAfter)
            {
                if (retryAfter < TimeSpan.Zero)
                    return TimeSpan.Zero;
                return retryAfter > MaxRetryAfter ? MaxRetryAfter : retryAfter;
            }

            return RetryBackoff[Math.Min(retry, RetryBackoff.Length - 1)];
        }

        /// <summary>The host a request URI goes to, resolving a relative one against the client's base address.</summary>
        internal string HostOf(string requestUri) =>
            Uri.TryCreate(requestUri, UriKind.Absolute, out Uri absolute)
                ? absolute.Host
                : _httpClient.BaseAddress?.Host ?? string.Empty;

        /// <summary>
        /// Names a download for an exception message by its file name and host only -- never the URL itself.
        /// </summary>
        /// <remarks>
        /// A download URL can carry a credential: PRIDE's reviewer-token route hands out hrefs whose query
        /// string IS the token. An exception message ends up in logs, CI output and pasted issues, so the URL
        /// never goes into one. The public route's URLs are harmless, but a rule that holds only for harmless
        /// URLs is not a rule.
        /// </remarks>
        internal static string DescribeDownload(string url, string fileName)
        {
            string host = Uri.TryCreate(url, UriKind.Absolute, out Uri uri) ? uri.Host : "an unparseable URL";
            return $"'{fileName}' from {host}";
        }

        /// <summary>
        /// One download attempt into <paramref name="partialPath"/>: a ranged request when the partial can be
        /// resumed, else a full one. Returns once the partial holds the complete file.
        /// </summary>
        /// <remarks>
        /// Bytes are appended only to a <c>206</c> whose <c>Content-Range</c> starts exactly where the partial
        /// ends, sent under <c>If-Range</c> so the server answers with the whole file instead if it has changed.
        /// With no validator there is no <c>If-Range</c>; the request is sent only when the MD5 checked afterwards
        /// will catch a changed file instead.
        /// Any other answer to a ranged request -- a <c>200</c>, a <c>206</c> starting elsewhere, a <c>416</c> --
        /// starts the file again from zero. The guarantee is that no truncated or spliced file is ever moved
        /// into place, not that no <c>.partial</c> is ever left behind.
        /// </remarks>
        internal async Task<bool> DownloadOnceAsync(string url, string described, string partialPath,
            string validatorPath, ResumableTransfer transfer, CancellationToken cancellationToken)
        {
            long resumeFrom = transfer.CanResume && File.Exists(partialPath) ? new FileInfo(partialPath).Length : 0;

            HttpResponseMessage response = await SendDownloadRequestAsync(url, described, resumeFrom, transfer.Validator, cancellationToken).ConfigureAwait(false);
            try
            {
                bool resumed = resumeFrom > 0
                    && response.StatusCode == HttpStatusCode.PartialContent
                    && response.Content.Headers.ContentRange?.From == resumeFrom;

                if (resumeFrom > 0 && !resumed
                    && response.StatusCode is HttpStatusCode.PartialContent or HttpStatusCode.RequestedRangeNotSatisfiable)
                {
                    // A 206 for some other range, or a 416, cannot be spliced onto what is on disk, and its body
                    // is not the whole file either. Ask again for the whole file.
                    response.Dispose();
                    resumeFrom = 0;
                    response = await SendDownloadRequestAsync(url, described, 0, null, cancellationToken).ConfigureAwait(false);
                }

                if (!response.IsSuccessStatusCode)
                    throw StatusFailure(response,
                        $"{_serviceName} download failed with status {(int)response.StatusCode} {response.ReasonPhrase} for {described}.");

                long? expectedLength;
                if (resumed)
                {
                    expectedLength = response.Content.Headers.ContentRange.Length;
                    transfer.Inspector?.ContinueAt(resumeFrom); // the start of the file is on disk, not in this stream
                }
                else
                {
                    // A fresh start: whatever is on disk is discarded, and the validator is the new response's.
                    expectedLength = response.Content.Headers.ContentLength;
                    transfer.Validator = ValidatorOf(response);
                    WriteValidator(validatorPath, transfer.Validator);
                    transfer.Inspector?.Reset();
                }

                using (var fileStream = new FileStream(partialPath, resumed ? FileMode.Append : FileMode.Create,
                           FileAccess.Write, FileShare.None))
                {
                    long writtenFrom = fileStream.Position; // the end of the partial when resumed, else 0
                    try
                    {
                        using (Stream httpStream = await response.Content.ReadAsStreamAsync(cancellationToken).ConfigureAwait(false))
                        {
                            // Not Stream.CopyToAsync: it would inherit the very absence of a read deadline that
                            // BodyStallTimeout exists to supply, which is how the body escaped every timeout here.
                            await CopyUntilStalledAsync(httpStream, fileStream, described, transfer.Inspector, cancellationToken).ConfigureAwait(false);
                        }
                    }
                    finally
                    {
                        transfer.BytesReceived += fileStream.Position - writtenFrom;
                    }
                }

                long actualLength = new FileInfo(partialPath).Length;
                if (expectedLength.HasValue && actualLength != expectedLength.Value)
                {
                    // The body ended cleanly but the file is not the size the server said it was. These bytes are
                    // not trusted for a resume either, so they go, and a retry starts from zero.
                    DeleteQuietly(partialPath);
                    DeleteQuietly(validatorPath);
                    transfer.Validator = null;
                    throw new HttpRequestException(
                        $"The {_serviceName} download of {described} came to {actualLength} bytes where the server announced {expectedLength.Value}.");
                }

                return true;
            }
            finally
            {
                response.Dispose();
            }
        }

        /// <summary>
        /// Sends the download GET. A positive <paramref name="resumeFrom"/> asks for the rest of the file only,
        /// and, given a <paramref name="validator"/>, only if the file still matches it.
        /// </summary>
        private Task<HttpResponseMessage> SendDownloadRequestAsync(string url, string described, long resumeFrom,
            RangeConditionHeaderValue validator, CancellationToken cancellationToken)
        {
            var request = new HttpRequestMessage(HttpMethod.Get, url);
            if (resumeFrom > 0)
            {
                request.Headers.Range = new RangeHeaderValue(resumeFrom, null);
                if (validator != null)
                    request.Headers.IfRange = validator;
            }

            return SendAsync(request, HttpCompletionOption.ResponseHeadersRead, cancellationToken, described);
        }

        /// <summary>
        /// The validator to resume a response's body under: its strong <c>ETag</c>, else its <c>Last-Modified</c>,
        /// else null. A weak ETag is not allowed in <c>If-Range</c>.
        /// </summary>
        /// <remarks>
        /// With no validator the bytes are resumed only when the whole file's MD5 is checked afterwards, not even
        /// within this call otherwise: a <c>Range</c> request with no <c>If-Range</c> could splice two versions of
        /// a file together, and only the MD5 would see it. PRIDE's public FTP host sends both validators. Its
        /// Globus HTTPS host sends neither, though it accepts <c>Range</c> (measured 2026-10-18), and neither does
        /// its reviewer-token route (measured 2026-09-29). There a download restarts from zero after a failure
        /// unless its MD5 is checked.
        /// </remarks>
        private static RangeConditionHeaderValue ValidatorOf(HttpResponseMessage response)
        {
            if (response.Headers.ETag is { IsWeak: false } etag)
                return new RangeConditionHeaderValue(etag);
            if (response.Content.Headers.LastModified is DateTimeOffset lastModified)
                return new RangeConditionHeaderValue(lastModified);
            return null;
        }

        /// <summary>
        /// Reads the validator a <c>.partial</c> was started under. An orphan -- a partial with no readable
        /// validator, or a validator with no partial -- is deleted, because it can never be resumed. With
        /// <paramref name="keepUnvalidated"/>, a partial with no validator is kept: the MD5 checked at the end of
        /// the download is what vouches for it then.
        /// </summary>
        internal static RangeConditionHeaderValue ReadValidator(string partialPath, string validatorPath, bool keepUnvalidated)
        {
            RangeConditionHeaderValue validator = null;
            try
            {
                if (File.Exists(partialPath) && File.Exists(validatorPath))
                    RangeConditionHeaderValue.TryParse(File.ReadAllText(validatorPath).Trim(), out validator);
            }
            catch (IOException) { }
            catch (UnauthorizedAccessException) { }

            if (validator == null)
            {
                if (!keepUnvalidated)
                    DeleteQuietly(partialPath);
                DeleteQuietly(validatorPath);
            }
            return validator;
        }

        /// <summary>Records <paramref name="validator"/> beside the partial, or removes a stale record when there is none.</summary>
        private static void WriteValidator(string validatorPath, RangeConditionHeaderValue validator)
        {
            if (validator == null)
                DeleteQuietly(validatorPath);
            else
                File.WriteAllText(validatorPath, validator.ToString());
        }

        /// <summary>
        /// Deletes a scratch file if it is there. Cleanup must never replace the exception that caused it: on
        /// Windows a locked file makes Delete throw, and ProteinDbRetriever.WriteResponseToFile guards the same way.
        /// </summary>
        internal static void DeleteQuietly(string path)
        {
            try
            {
                if (File.Exists(path))
                    File.Delete(path);
            }
            catch (IOException) { }
            catch (UnauthorizedAccessException) { }
        }

        /// <summary>
        /// Copies <paramref name="source"/> to <paramref name="destination"/>, giving up if no bytes arrive
        /// for <see cref="BodyStallTimeout"/>. Each read gets a FRESH window, so a transfer that keeps
        /// delivering data runs as long as it needs to and only an actual stoppage ends it.
        /// </summary>
        /// <remarks>
        /// The stall window is linked to <paramref name="cancellationToken"/> so a caller's cancellation is
        /// still honoured immediately, but the two are told apart afterwards: only a window that fired on
        /// its own becomes an <see cref="HttpRequestException"/>. A caller who cancels gets the
        /// <see cref="OperationCanceledException"/> they asked for, because reporting that as a transport
        /// failure would make <c>ExternalServiceTestHelper</c> skip a test that was deliberately cancelled.
        /// <para>
        /// Internal, not private, because it is the one part of a download that is not HTTP's: the MassIVE client
        /// copies an FTP data connection through it, whose read failures arrive as an <see cref="IOException"/> too.
        /// </para>
        /// </remarks>
        internal async Task CopyUntilStalledAsync(Stream source, Stream destination, string described,
            IDownloadInspector inspector, CancellationToken cancellationToken)
        {
            byte[] buffer = new byte[81920];
            while (true)
            {
                int read;
                using (var stallWindow = new CancellationTokenSource(BodyStallTimeout))
                using (var linked = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken, stallWindow.Token))
                {
                    try
                    {
                        read = await source.ReadAsync(buffer.AsMemory(), linked.Token).ConfigureAwait(false);
                    }
                    catch (OperationCanceledException e) when (stallWindow.IsCancellationRequested
                                                               && !cancellationToken.IsCancellationRequested)
                    {
                        throw new HttpRequestException(
                            $"The {_serviceName} response body for {described} delivered nothing for {BodyStallTimeout}.", e);
                    }
                    catch (IOException e)
                    {
                        // A connection the server drops mid-body arrives here as an IOException (HttpIOException on
                        // .NET 10, "The response ended prematurely"), not as the HttpRequestException the
                        // contract promises -- so a caller catching "try again later" missed the phase where a
                        // large download actually breaks, and a live test reddened instead of skipping. Only
                        // the READ is wrapped: a failure writing the local file is the caller's disk, not an
                        // outage, and stays an IOException.
                        throw new HttpRequestException(
                            $"The {_serviceName} response body for {described} ended before the download completed.", e);
                    }
                }

                if (read == 0)
                    return;

                inspector?.Append(buffer.AsSpan(0, read));
                if (inspector?.StopReason is string stopReason)
                {
                    // The inspector has seen enough to know this attempt can only fail: stop paying for the rest of it.
                    // A damaged 2.3 GB copy cost 10-40 minutes when it was downloaded to the end (PXReprise 019).
                    throw new ContentMismatchException($"The {_serviceName} download of {described} {stopReason}.");
                }
                await destination.WriteAsync(buffer.AsMemory(0, read), cancellationToken).ConfigureAwait(false);
            }
        }

        /// <summary>
        /// Sends a GET for one byte range of a file and returns the response, its body not yet read, once it is a
        /// <c>206</c> for exactly that range: bytes <paramref name="from"/> to <paramref name="to"/>, or, with
        /// <paramref name="from"/> null, the last <paramref name="to"/> bytes (the whole file when it is shorter).
        /// One attempt, no retry: the caller's attempt, inside <see cref="WithRetryAsync{T}"/>, is the unit retried.
        /// </summary>
        /// <remarks>
        /// <para>
        /// A non-success status throws <see cref="StatusFailure"/>, worded as a download's ("failed with status NNN").
        /// </para>
        /// <para>
        /// A <c>200</c> is never read: it is the whole file, which for a large one is the very transfer a range read
        /// exists to avoid. Under <paramref name="validator"/> (sent as <c>If-Range</c>) a <c>200</c> means the file
        /// changed since the validator was taken (PRIDE's FTP host, measured 2026-11-18: a wrong <c>If-Range</c>
        /// answers <c>200</c> with all 25 MB of the file); without one, that the server ignores ranges. Neither is an
        /// outage, so it is an <see cref="MzLibException"/>, and so is a <c>Content-Range</c> total other than
        /// <paramref name="expectedTotal"/>: the file changed under a reader with no validator.
        /// </para>
        /// <para>
        /// A <c>206</c> for some other range is a <see cref="ContentMismatchException"/>, retried as damaged bytes are.
        /// </para>
        /// </remarks>
        /// <param name="url">The file's URL. It never enters a message.</param>
        /// <param name="described">How messages name the file (<see cref="DescribeDownload"/>).</param>
        /// <param name="from">The first byte wanted, or null for the last <paramref name="to"/> bytes.</param>
        /// <param name="to">The last byte wanted (inclusive), or with <paramref name="from"/> null the suffix's length.</param>
        /// <param name="validator">Sent as <c>If-Range</c> when not null.</param>
        /// <param name="expectedTotal">The file's length as an earlier read reported it, or null on the first read.</param>
        /// <param name="cancellationToken">The caller's token.</param>
        internal async Task<HttpResponseMessage> GetRangeAsync(string url, string described, long? from, long to,
            RangeConditionHeaderValue validator, long? expectedTotal, CancellationToken cancellationToken)
        {
            var request = new HttpRequestMessage(HttpMethod.Get, url);
            request.Headers.Range = new RangeHeaderValue(from, to);
            if (validator != null)
                request.Headers.IfRange = validator;

            HttpResponseMessage response = await SendAsync(request, HttpCompletionOption.ResponseHeadersRead, cancellationToken, described).ConfigureAwait(false);
            try
            {
                if (!response.IsSuccessStatusCode)
                    throw StatusFailure(response,
                        $"{_serviceName} download failed with status {(int)response.StatusCode} {response.ReasonPhrase} for {described}.");

                if (response.StatusCode != HttpStatusCode.PartialContent)
                    throw new MzLibUtil.MzLibException(validator != null
                        ? $"{_serviceName} answered a request for part of {described} with the whole file: the file has changed since its first part was read, or the server ignores byte ranges."
                        : $"{_serviceName} answered a request for part of {described} with the whole file: the server ignores byte ranges.");

                ContentRangeHeaderValue served = response.Content.Headers.ContentRange;
                if (served is not { HasRange: true, HasLength: true } || served.Unit != "bytes")
                    throw new ContentMismatchException($"The {_serviceName} answer for part of {described} did not say which bytes it holds.");

                if (expectedTotal.HasValue && served.Length != expectedTotal.Value)
                    throw new MzLibUtil.MzLibException(
                        $"{described} changed on the server while it was being read: it is now {served.Length} bytes, not {expectedTotal.Value}.");

                long wantedFrom = from ?? Math.Max(0, served.Length.Value - to);
                long wantedTo = from.HasValue ? to : served.Length.Value - 1;
                if (served.From != wantedFrom || served.To != wantedTo)
                    throw new ContentMismatchException(
                        $"The {_serviceName} answer for bytes {wantedFrom}-{wantedTo} of {described} held bytes {served.From}-{served.To} instead.");

                return response;
            }
            catch
            {
                response.Dispose();
                throw;
            }
        }

        /// <summary>
        /// <see cref="GetRangeAsync"/>, with the body read into memory under <see cref="BodyStallTimeout"/>, and checked to be
        /// as long as its <c>Content-Range</c> says. One attempt, no retry.
        /// </summary>
        internal async Task<RangeBytes> ReadRangeOnceAsync(string url, string described, long? from, long to,
            RangeConditionHeaderValue validator, long? expectedTotal, CancellationToken cancellationToken)
        {
            using HttpResponseMessage response = await GetRangeAsync(url, described, from, to, validator, expectedTotal, cancellationToken).ConfigureAwait(false);
            ContentRangeHeaderValue served = response.Content.Headers.ContentRange;
            long length = served.To.Value - served.From.Value + 1;
            if (length > int.MaxValue)
                throw new NotSupportedException($"A {length}-byte part of {described} is too large to read into memory.");

            using var bytes = new MemoryStream((int)length);
            using (Stream body = await response.Content.ReadAsStreamAsync(cancellationToken).ConfigureAwait(false))
                await CopyUntilStalledAsync(body, bytes, described, new LengthLimit(length), cancellationToken).ConfigureAwait(false);

            if (bytes.Length != length)
                throw new HttpRequestException(
                    $"The {_serviceName} answer for part of {described} came to {bytes.Length} bytes where the server announced {length}.");

            return new RangeBytes(bytes.ToArray(), served.From.Value, served.Length.Value, ValidatorOf(response));
        }

        /// <summary>One byte range of a file, read whole.</summary>
        /// <param name="Bytes">The bytes.</param>
        /// <param name="From">The offset of the first of them in the file.</param>
        /// <param name="Total">The file's length, from <c>Content-Range</c>.</param>
        /// <param name="Validator">
        /// The validator later reads of the same file are sent under (strong <c>ETag</c>, else <c>Last-Modified</c>), or
        /// null when the server sent neither.
        /// </param>
        internal sealed record RangeBytes(byte[] Bytes, long From, long Total, RangeConditionHeaderValue Validator);

        /// <summary>Stops a read that brings more bytes than its range holds, instead of buffering them all.</summary>
        private sealed class LengthLimit(long limit) : IDownloadInspector
        {
            private long _seen;

            public void Reset() => _seen = 0;

            public void ContinueAt(long offset) => _seen = offset;

            public void Append(ReadOnlySpan<byte> data) => _seen += data.Length;

            public string StopReason => _seen > limit ? $"brought more than the {limit} bytes asked for" : null;
        }
    }

    /// <summary>
    /// Watches a download's bytes as they arrive, for a client that checks them against something of its own
    /// (PRIDE: its checksum list's MD5, or EBI's error text inside the file).
    /// </summary>
    internal interface IDownloadInspector
    {
        /// <summary>Starts over, for a file that is about to arrive from its first byte.</summary>
        void Reset();

        /// <summary>Continues a resumed transfer whose bytes resume at <paramref name="offset"/> in the file.</summary>
        void ContinueAt(long offset);

        /// <summary>Takes in the next bytes of the file, in order.</summary>
        void Append(ReadOnlySpan<byte> data);

        /// <summary>
        /// Null while the attempt may still succeed; once the bytes seen show it cannot, the end of a sentence saying
        /// why, naming no path or URL. The attempt is then stopped as a <see cref="ContentMismatchException"/>.
        /// </summary>
        string StopReason { get; }
    }

    /// <summary>What one download carries from a failed attempt into the next.</summary>
    internal sealed class ResumableTransfer
    {
        /// <summary>The validator the partial bytes were served under, or null when there is none.</summary>
        public RangeConditionHeaderValue Validator { get; set; }

        /// <summary>
        /// Whether the partial bytes may be resumed with no validator: only when the whole file's MD5 is checked
        /// at the end, which a splice of two versions of the file fails.
        /// </summary>
        public bool ResumableWithoutValidator { get; init; }

        /// <summary>Whether the partial bytes may be resumed at all.</summary>
        public bool CanResume => Validator != null || ResumableWithoutValidator;

        /// <summary>The check run over the bytes as they arrive, or null when there is nothing to check them against.</summary>
        public IDownloadInspector Inspector { get; init; }

        /// <summary>The bytes the attempts so far have written to the partial, counted as they arrive, kept or not.</summary>
        public long BytesReceived { get; set; }
    }

    /// <summary>
    /// A downloaded file is not the attempt's failure but its content's: wrong size, wrong MD5, or a server's error
    /// text inside it. It has no status, so <see cref="DownloadEngine.WithRetryAsync{T}"/> treats it as transient and
    /// downloads the file again from zero; once the retries are spent the client turns it into an
    /// <c>MzLibException</c>.
    /// </summary>
    internal sealed class ContentMismatchException(string message) : HttpRequestException(message);
}
