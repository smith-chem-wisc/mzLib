using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Security.Cryptography;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A client for the PRIDE Archive REST API (https://www.ebi.ac.uk/pride/ws/archive/v3/) that makes
    /// EBI PRIDE proteomics dataset data available to mzLib. Given a project accession (e.g. "PXD012345")
    /// it returns the project's metadata, the complete manifest of its files, and the files' bytes; it
    /// also fetches individual spectra by USI through the PROXI standard.
    /// </summary>
    /// <remarks>
    /// Two failure kinds are deliberately given different exception types, because callers — and
    /// <c>ExternalServiceTestHelper</c>, which skips a live test on <see cref="HttpRequestException"/> —
    /// must be able to tell them apart. The service being unreachable, timing out or answering 5xx is an
    /// <see cref="HttpRequestException"/>. PRIDE answering successfully but with something the contract
    /// forbids — no such project, an empty body, a payload missing its accession — is an
    /// <see cref="MzLibException"/>, so it fails loudly instead of being mistaken for an outage.
    /// <para>
    /// Holds one reusable <see cref="HttpClient"/> and is <see cref="IDisposable"/>. Use one instance
    /// per unit of work and dispose it. A constructor overload accepts an <see cref="HttpClient"/> to
    /// support testing and custom configuration.
    /// </para>
    /// </remarks>
    public sealed class PrideArchiveClient : IDisposable
    {
        /// <summary>The base address of the PRIDE Archive REST API (v3).</summary>
        public const string DefaultBaseAddress = "https://www.ebi.ac.uk/pride/ws/archive/v3/";

        /// <summary>
        /// The base address of PRIDE's PROXI spectrum API. PROXI (the PSI standard for retrieving a spectrum by
        /// USI) lives under a DIFFERENT path root than the v3 archive API (<see cref="DefaultBaseAddress"/>), so
        /// PROXI requests are issued as absolute URIs and do not resolve against the client's archive
        /// <see cref="HttpClient.BaseAddress"/>.
        /// </summary>
        public const string DefaultProxiBaseAddress = "https://www.ebi.ac.uk/pride/proxi/archive/v0.1/";

        // Explicit JSON nulls must not clobber the non-null string defaults on the DTOs.
        private static readonly JsonSerializerSettings JsonSettings = new() { NullValueHandling = NullValueHandling.Ignore };

        private readonly HttpClient _httpClient;
        private readonly bool _ownsHttpClient;
        private bool _disposed;

        /// <summary>
        /// A hard upper bound on the number of pages fetched, guarding against a misbehaving server that
        /// ignores the paging parameters. No real PRIDE project approaches this; exceeding it throws.
        /// </summary>
        public int MaxPages { get; init; } = 10000;

        /// <summary>
        /// The longest keyword <see cref="SearchProjectsAsync"/> will send. PRIDE answers a very long
        /// keyword with HTTP 500 rather than a 400 or a 414 (observed at 2000 characters on
        /// 2026-08-21; 500 characters still answered 200), and a 500 is the one signature that cannot
        /// be told apart from an outage — <c>ExternalServiceTestHelper</c> reads it as "the service is
        /// down" and SKIPS. Refusing the request here turns a caller's bug into an
        /// <see cref="ArgumentException"/> at the call site instead. The exact server threshold is not
        /// published; this sits well inside the range observed to work.
        /// </summary>
        public const int MaxKeywordLength = 1000;

        /// <summary>
        /// How long <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/> will wait for the NEXT bytes of a response body before
        /// abandoning the transfer. Each read gets a fresh window, so this bounds silence, not duration.
        /// </summary>
        /// <remarks>
        /// An INACTIVITY deadline, deliberately not a total-duration one. PRIDE serves raw and peak files
        /// measured in gigabytes and a slow link may legitimately need a very long time to finish them, so
        /// the question worth asking is whether data is still arriving — not how long it has been arriving
        /// for. A total cap would abort exactly the healthy downloads that need the most time.
        /// <para>
        /// <see cref="HttpClient.Timeout"/> does not cover this. The download reads with
        /// <see cref="HttpCompletionOption.ResponseHeadersRead"/>, so the client timeout is spent once the
        /// headers arrive and the body that follows had no deadline at all — a stalled transfer stayed open
        /// until the socket gave up.
        /// </para>
        /// <para>
        /// The failure this prevents was observed in this repository, in the sibling client rather than
        /// here: on 2026-08-21 a stalled UniProt body occupied 12m43s of a 20-minute CI job before the job
        /// was cancelled, reddening the run on every open PR. <c>ProteinDbRetriever.BodyStallTimeout</c>
        /// (#1189) closed it there; this closes the same hole here, and the two minutes is that fix's
        /// measured value, not a fresh guess — a merely SLOW transfer is unaffected because every byte
        /// restarts the clock, and two minutes of complete silence on a response body is already a dead
        /// connection.
        /// </para>
        /// <para>
        /// A stall is reported as <see cref="HttpRequestException"/>, the transport-failure type: silence
        /// from EBI is an outage, not a contract break, and a live test should skip on it rather than fail.
        /// Cancellation by the caller stays an <see cref="OperationCanceledException"/> and is not converted.
        /// </para>
        /// <para>
        /// Settable through an object initializer, as <see cref="MaxPages"/> is, so a test can prove a stall
        /// is detected without waiting minutes for it.
        /// </para>
        /// </remarks>
        public TimeSpan BodyStallTimeout { get; init; } = TimeSpan.FromMinutes(2);

        /// <summary>
        /// The wait before each retry of a transient failure: 5 s, 20 s, then 60 s, so a failure that persists
        /// gives up after about 85 s of backoff. See <see cref="WithRetryAsync{T}"/> for what counts as transient.
        /// </summary>
        /// <remarks>
        /// Deliberately internal. mzLib has no public retry policy anywhere, and nobody has asked for one;
        /// a knob can be published later, but a published one cannot be taken back.
        /// </remarks>
        private static readonly TimeSpan[] RetryBackoff =
            { TimeSpan.FromSeconds(5), TimeSpan.FromSeconds(20), TimeSpan.FromSeconds(60) };

        /// <summary>The longest <c>Retry-After</c> a 429 or 503 is allowed to make a retry wait.</summary>
        private static readonly TimeSpan MaxRetryAfter = TimeSpan.FromSeconds(60);

        /// <summary>How many times a transient failure is retried after the first attempt. Tests set 0 to pin a single attempt.</summary>
        internal int MaxRetries { get; init; } = RetryBackoff.Length;

        /// <summary>Waits out a backoff. Tests replace it to record the schedule instead of sleeping through it.</summary>
        internal Func<TimeSpan, CancellationToken, Task> RetryDelay { get; init; } = Task.Delay;

        /// <summary>Creates a client with its own <see cref="HttpClient"/> pointed at the PRIDE Archive API.</summary>
        public PrideArchiveClient()
            : this(new HttpClient { BaseAddress = new Uri(DefaultBaseAddress), Timeout = TimeSpan.FromSeconds(100) }, ownsHttpClient: true)
        {
        }

        /// <summary>
        /// Creates a client over a supplied <see cref="HttpClient"/> (for testing or custom configuration).
        /// The caller retains ownership: the supplied client is NOT disposed by <see cref="Dispose"/>.
        /// If the client has no <see cref="HttpClient.BaseAddress"/>, the PRIDE Archive base address is set.
        /// </summary>
        /// <param name="httpClient">The HTTP client to use. Must not be null.</param>
        public PrideArchiveClient(HttpClient httpClient)
            : this(httpClient, ownsHttpClient: false)
        {
        }

        private PrideArchiveClient(HttpClient httpClient, bool ownsHttpClient)
        {
            _httpClient = httpClient ?? throw new ArgumentNullException(nameof(httpClient));
            _ownsHttpClient = ownsHttpClient;
            if (_httpClient.BaseAddress == null)
                _httpClient.BaseAddress = new Uri(DefaultBaseAddress);
        }

        /// <summary>
        /// Returns the manifest of files belonging to a PRIDE Archive project. All pages are fetched
        /// and concatenated; paging is an implementation detail hidden from the caller. The manifest is
        /// complete at the default page size and for any size at or below the server cap, and complete
        /// above the cap as long as the response reports <c>total_records</c> (see
        /// <paramref name="pageSize"/> for the one case that can still return a partial manifest).
        /// </summary>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="pageSize">
        /// Files requested per page (default 100). PRIDE silently caps this server-side and then pages
        /// by the capped size, so a value above the cap means more pages rather than fewer files
        /// <em>provided the response reports <c>total_records</c></em>, which PRIDE does on every
        /// observed response. If that header is ever absent or unparseable there is nothing left to
        /// page against, and termination falls back to the first short page — under which a requested
        /// size above the server cap can still return a partial manifest. Leave this at the default
        /// unless you have a reason not to; the full manifest is returned either way at or below the cap.
        /// <para>
        /// A <c>total_records</c> that misreports the count in either direction is tolerated: an
        /// overstated one stops on the empty page past the end, and an understated one is caught by
        /// fetching one further page whenever the fetch would otherwise end on a page that could be
        /// full — a first page holding everything that was asked for, or a later page as large as the
        /// largest the server has served. That costs one extra request per fetch of exactly that shape;
        /// a fetch ending on a page that is visibly short pays nothing. The one case it cannot catch is
        /// an understated total on a single page requested ABOVE the server cap, where the capped page
        /// that comes back looks short but may not be — the same above-the-cap caveat as the paragraph
        /// above.
        /// </para>
        /// </param>
        /// <param name="cancellationToken">Cancels the (possibly multi-page) fetch.</param>
        /// <returns>
        /// The project's files. Empty if the project has no files or the accession is unknown (PRIDE
        /// returns an empty result for an unknown accession). Never null.
        /// </returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The page size is not positive.</exception>
        /// <exception cref="HttpRequestException">The API returned a non-success status code or did not answer. A transient failure of one page is retried, that page alone, as <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/> describes before this is thrown.</exception>
        /// <exception cref="MzLibException">
        /// PRIDE answered successfully but served a page identical to its predecessor while
        /// <c>total_records</c> reported more remained — a broken contract rather than an outage, so it
        /// is not mistaken for one (see the class remarks) and is not retried as network trouble.
        /// </exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<List<PrideArchiveFile>> GetProjectFilesAsync(string accession, int pageSize = 100,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A PRIDE project accession is required.", nameof(accession));
            if (pageSize <= 0)
                throw new ArgumentOutOfRangeException(nameof(pageSize), "Page size must be positive.");

            return await GetAllPagesAsync<PrideArchiveFile>(
                page => $"projects/{Uri.EscapeDataString(accession)}/files?pageSize={pageSize}&page={page}",
                pageSize,
                $"accession '{accession}'",
                cancellationToken).ConfigureAwait(false);
        }

        /// <summary>
        /// Returns the COMPLETE list of a project's files by walking its FTP directory tree — for the
        /// cases where <see cref="GetProjectFilesAsync"/> is not enough. PRIDE's REST manifest is
        /// knowingly incomplete (for PXD000001 it lists 8 files while the FTP tree holds 13 — it omits
        /// five, including the two largest), so a caller that must know everything a project contains —
        /// or its true size — reads the tree instead of trusting the manifest.
        /// </summary>
        /// <remarks>
        /// The tree lives at <c>https://{PrideFtpHost}/pride/data/archive/{yyyy}/{MM}/{accession}/</c>,
        /// where the year and month are the project's <see cref="PrideProject.PublicationDate"/>. The FTP
        /// host also serves that path over HTTPS — the same fact <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/> relies
        /// on — so the whole walk goes over this client's reused <see cref="HttpClient"/>. Subdirectories
        /// are followed; each returned file's <see cref="PrideFtpFile.RelativePath"/> is relative to the
        /// project root. Sizes are PRIDE's rounded index sizes — see <see cref="PrideFtpFile.ApproximateSizeBytes"/>.
        /// </remarks>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="cancellationToken">Cancels the (multi-request) walk.</param>
        /// <returns>Every file under the project's FTP root, subdirectories included. Never null.</returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="MzLibException">No project has that accession, or it carries no publication date to locate its FTP directory.</exception>
        /// <exception cref="HttpRequestException">The project or a directory could not be fetched (non-success status), or the project's root directory listed no entries at all, which is EBI serving a bad listing rather than an empty project. Transient failures, an empty root included, are retried as <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/> describes before this is thrown.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<List<PrideFtpFile>> GetProjectFilesFromFtpAsync(string accession,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A PRIDE project accession is required.", nameof(accession));

            // GetProjectAsync gives the publication date AND the right failure contract: an unknown
            // accession is an MzLibException (not an HttpRequestException that a live test would skip).
            PrideProject project = await GetProjectAsync(accession, cancellationToken).ConfigureAwait(false);
            if (project.PublicationDate == default)
                throw new MzLibException(
                    $"PRIDE project '{accession}' has no publication date, so its FTP directory cannot be located.");

            string rootUrl = string.Format(CultureInfo.InvariantCulture,
                "https://{0}/pride/data/archive/{1:yyyy}/{1:MM}/{2}/",
                PrideArchiveExtensions.PrideFtpHost, project.PublicationDate, accession);

            var files = new List<PrideFtpFile>();
            await CollectFtpFilesAsync(accession, new Uri(rootUrl), relativePrefix: "", depth: 0, files, cancellationToken)
                .ConfigureAwait(false);
            return files;
        }

        /// <summary>
        /// A hard cap on how deep the FTP walk descends. A symlink cycle is normally pruned by the
        /// visited-URL set, so this is only a last-resort backstop for a genuinely (and implausibly)
        /// deep tree; no real PRIDE project nests anywhere near this deep, and exceeding it is treated
        /// as a broken listing (<see cref="MzLibException"/>) rather than left to a process-killing
        /// <see cref="StackOverflowException"/>.
        /// </summary>
        private const int MaxFtpDirectoryDepth = 64;

        /// <summary>
        /// Lists one FTP directory (over HTTPS), appends its files to <paramref name="files"/>, and
        /// recurses into each child subdirectory. <paramref name="relativePrefix"/> is the path from the
        /// project root down to this directory, so nested files carry a full relative path;
        /// <paramref name="depth"/> bounds the recursion.
        /// </summary>
        private async Task CollectFtpFilesAsync(string accession, Uri directoryUri, string relativePrefix, int depth,
            List<PrideFtpFile> files, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();

            // The traversal guard below only ever appends a single non-".." segment, so the URL strictly
            // deepens and can never revisit an earlier directory — no visited-set is needed. A cyclic
            // listing (a self-linking symlink) therefore grows the URL without bound and is caught here.
            // A broken listing is a contract violation (MzLibException) — NOT an HttpRequestException,
            // which ExternalServiceTestHelper would misread as a service outage.
            if (depth > MaxFtpDirectoryDepth)
                throw new MzLibException(
                    $"PRIDE FTP directory nesting exceeded {MaxFtpDirectoryDepth} levels at '{directoryUri}'; the listing may be cyclic.");

            List<(string RawHref, string SizeText)> rows = await WithRetryAsync(async () =>
            {
                using HttpResponseMessage response = await GetAsync(directoryUri.AbsoluteUri, HttpCompletionOption.ResponseContentRead, cancellationToken).ConfigureAwait(false);
                if (!response.IsSuccessStatusCode)
                    throw StatusFailure(response,
                        $"PRIDE FTP directory listing failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{directoryUri}'.");

                string html = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);
                List<(string RawHref, string SizeText)> listed = ParseAutoIndexRows(html).ToList();

                // The project resolved, so its root directory exists and holds at least the files PRIDE was
                // given. A root that answers 200 with no entries is EBI serving a bad listing, and returning []
                // for it made PXD058248 (46 files) look like "no such project" to pyMzLib, which does not retry
                // that. So it is a transport failure, and retried. A SUBDIRECTORY can genuinely be empty
                // (PXD001174's .wiff directory is), so only the root is held to this.
                if (depth == 0 && listed.Count == 0)
                    throw new HttpRequestException(
                        $"PRIDE FTP listing of project {accession}'s root directory came back empty although the project exists.");

                return listed;
            }, directoryUri.Host, cancellationToken).ConfigureAwait(false);

            foreach ((string rawHref, string sizeText) in rows)
            {
                // rawHref is the link's attribute value, still URL-encoded. Resolve HTML entities, then
                // compose the child URL with Uri joining (which escapes and handles the base's trailing
                // slash safely), and decode only for the human-readable name — so "my%20file.raw"
                // downloads from an encoded URL but surfaces as "my file.raw".
                string href = WebUtility.HtmlDecode(rawHref);
                string name = Uri.UnescapeDataString(href);
                bool isDirectory = href.EndsWith("/", StringComparison.Ordinal);
                string segment = isDirectory ? name.TrimEnd('/') : name;

                // Refuse a "." / ".." traversal (which the Uri join would normalise OUT of the project
                // subtree) and any decoded segment that gained a path separator (an encoded "%2F"),
                // before it is used to build either the URL or the relative path.
                if (segment.Length == 0 || segment == "." || segment == ".." || segment.Contains('/'))
                    continue;

                Uri childUri = new(directoryUri, href);

                if (isDirectory)
                    await CollectFtpFilesAsync(accession, childUri, relativePrefix + name, depth + 1, files, cancellationToken)
                        .ConfigureAwait(false);
                else
                    files.Add(new PrideFtpFile(relativePrefix + name, ParseAutoIndexSize(sizeText), childUri.AbsoluteUri));
            }
        }

        // The PRIDE FTP host serves a standard Apache autoindex table over HTTPS: one <tr> per entry,
        // whose <a href> is the (relative) name and whose right-aligned cells are the last-modified date
        // and the size. The sort headers link to "?C=...", and "Parent Directory" links to an absolute
        // "/..." path; both are skipped so only real entries survive.
        private static readonly Regex RowRegex = new(@"<tr[^>]*>(.*?)</tr>", RegexOptions.Singleline | RegexOptions.IgnoreCase);
        private static readonly Regex HrefRegex = new("<a\\s+href=\"([^\"]+)\"", RegexOptions.IgnoreCase);
        private static readonly Regex RightCellRegex = new("<td[^>]*align=\"right\"[^>]*>(.*?)</td>", RegexOptions.Singleline | RegexOptions.IgnoreCase);
        // A size cell is "-" (a directory) or a number optionally suffixed K/M/G/T — never a date.
        private static readonly Regex SizeRegex = new(@"^(-|\d+(\.\d+)?[KMGT]?)$", RegexOptions.IgnoreCase);

        /// <summary>Yields the raw href and the size text for each real entry in an Apache autoindex table.</summary>
        private static IEnumerable<(string RawHref, string SizeText)> ParseAutoIndexRows(string html)
        {
            foreach (Match row in RowRegex.Matches(html))
            {
                Match href = HrefRegex.Match(row.Groups[1].Value);
                if (!href.Success)
                    continue;

                string raw = href.Groups[1].Value;
                // Skip the column-sort links ("?C=N;O=D"), the absolute-path "Parent Directory", and
                // anchors. Test the leading char on the decoded form; the RAW href is what the caller joins.
                string decoded = WebUtility.HtmlDecode(raw);
                if (decoded.Length == 0 || decoded[0] == '/' || decoded[0] == '?' || decoded[0] == '#')
                    continue;

                // Pick the FIRST right-aligned cell whose text parses as a size, NOT blindly the last: a
                // stock Apache template also right-aligns the (empty) Description column, and the date cell
                // is right-aligned too — neither of which looks like a size.
                string sizeText = string.Empty;
                foreach (Match cell in RightCellRegex.Matches(row.Groups[1].Value))
                {
                    string text = WebUtility.HtmlDecode(cell.Groups[1].Value).Trim();
                    if (SizeRegex.IsMatch(text))
                    {
                        sizeText = text;
                        break;
                    }
                }

                yield return (raw, sizeText);
            }
        }

        /// <summary>
        /// Parses one Apache autoindex size cell ("1.6K", "20M", "9.3M", "1.2G", a bare byte count, or
        /// "-" for a directory) into an APPROXIMATE byte count. The index rounds to ~3 significant
        /// figures, so this is not exact — see <see cref="PrideFtpFile.ApproximateSizeBytes"/>.
        /// </summary>
        private static long ParseAutoIndexSize(string sizeText)
        {
            sizeText = sizeText.Trim();
            if (sizeText.Length == 0 || sizeText == "-")
                return 0;

            char last = char.ToUpperInvariant(sizeText[sizeText.Length - 1]);
            bool hasSuffix = last is 'K' or 'M' or 'G' or 'T';
            double multiplier = last switch
            {
                'K' => 1024d,
                'M' => 1024d * 1024,
                'G' => 1024d * 1024 * 1024,
                'T' => 1024d * 1024 * 1024 * 1024,
                _ => 1d,
            };
            string number = hasSuffix ? sizeText.Substring(0, sizeText.Length - 1) : sizeText;
            return double.TryParse(number, NumberStyles.Float, CultureInfo.InvariantCulture, out double value)
                ? (long)Math.Round(value * multiplier)
                : 0;
        }

        /// <summary>
        /// Returns the metadata describing a PRIDE Archive project — title, protocols, instruments,
        /// species, publications and submitters — letting a caller judge and cite a dataset before
        /// downloading any of it. Use <see cref="TryGetProjectAsync"/> instead when the accession comes
        /// from user input and may simply not exist.
        /// </summary>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns>The project's metadata. Never null.</returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="MzLibException">
        /// No project has that accession — unlike <see cref="GetProjectFilesAsync"/>, which answers an
        /// unknown accession with an empty manifest, this endpoint answers with 404. This is deliberately
        /// NOT an <see cref="HttpRequestException"/>: that type means "the service is unavailable" and is
        /// converted to a skipped test by <c>ExternalServiceTestHelper</c>, which would let a withdrawn
        /// accession pass unnoticed instead of failing.
        /// </exception>
        /// <exception cref="HttpRequestException">The API was unreachable or returned a non-success status other than 404.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<PrideProject> GetProjectAsync(string accession, CancellationToken cancellationToken = default)
        {
            (bool found, PrideProject project) = await TryGetProjectAsync(accession, cancellationToken).ConfigureAwait(false);

            if (!found)
                throw new MzLibException($"PRIDE Archive has no project with accession '{accession}'.");

            return project;
        }

        /// <summary>
        /// Attempts to fetch a PRIDE Archive project's metadata, reporting a non-existent accession as a
        /// value rather than an exception — the cheap way to validate an accession a user typed.
        /// </summary>
        /// <remarks>
        /// <c>Found</c> is false for one reason only: PRIDE answered, and no project has that accession
        /// (HTTP 404). Every other failure — the service being unreachable, timing out, rate-limiting, or
        /// returning 5xx — still throws, because collapsing an outage into "no such project" would send a
        /// caller hunting for a typo in a perfectly good accession. This is also what lets a live test
        /// distinguish "PRIDE is down, skip" from "the contract broke, fail".
        /// </remarks>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns>
        /// <c>Found</c> and the project when the accession exists; <c>false</c> and null when PRIDE reports
        /// no such project. The project is never null when <c>Found</c> is true.
        /// </returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="HttpRequestException">The API was unreachable or returned a non-success status other than 404.</exception>
        /// <exception cref="MzLibException">PRIDE answered successfully but the payload was empty or carried no accession — a broken contract rather than an absence.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<(bool Found, PrideProject Project)> TryGetProjectAsync(string accession,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A PRIDE project accession is required.", nameof(accession));

            cancellationToken.ThrowIfCancellationRequested();

            // This endpoint shares the v3 BaseAddress, so a relative URI resolves correctly (unlike PROXI,
            // which sits under a different path root and needs an absolute URI).
            string requestUri = $"projects/{Uri.EscapeDataString(accession)}";
            using HttpResponseMessage response = await WithRetryAsync(async () =>
            {
                HttpResponseMessage attempt = await GetAsync(requestUri, HttpCompletionOption.ResponseContentRead, cancellationToken).ConfigureAwait(false);

                // 404 is the ONE expected "no" and is reported as a value, so it is let through here and never
                // retried. It is checked before the general status guard so that every other failure still throws.
                if (attempt.IsSuccessStatusCode || attempt.StatusCode == HttpStatusCode.NotFound)
                    return attempt;

                using (attempt)
                    throw StatusFailure(attempt,
                        $"PRIDE Archive request failed with status {(int)attempt.StatusCode} {attempt.ReasonPhrase} for '{requestUri}'.");
            }, HostOf(requestUri), cancellationToken).ConfigureAwait(false);

            if (response.StatusCode == HttpStatusCode.NotFound)
                return (false, null);

            string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

            // Unlike the manifest and PROXI endpoints, this one returns a bare object rather than an array,
            // so there is nothing to unwrap. A 200 carrying no project is a broken contract, not an absence,
            // so it throws MzLibException rather than HttpRequestException — the latter would be read as an
            // outage and silently skip a live test.
            PrideProject project = JsonConvert.DeserializeObject<PrideProject>(content, JsonSettings);
            if (project == null)
                throw new MzLibException($"PRIDE Archive returned an empty body for accession '{accession}'.");

            // An empty JSON object deserializes to a fully-defaulted instance, which would otherwise be
            // handed back as a successful result. The accession is always present on a real project, so its
            // absence means the payload is not one.
            if (string.IsNullOrEmpty(project.Accession))
                throw new MzLibException(
                    $"PRIDE Archive returned a payload with no accession for '{accession}'; the response is not a project.");

            // NullValueHandling.Ignore suppresses null *values*, not null *elements*: PRIDE sending
            // "instruments": [null] would otherwise put a null inside a collection the DTO documents as
            // never-null, and NRE at the caller's first dereference.
            RemoveNullElements(project);

            return (true, project);
        }

        /// <summary>
        /// Returns the MD5 and exact size PRIDE recorded for each file deposited in a project, keyed by bare
        /// file name — what a caller needs to check a download, or to tell whether a file already on disk is
        /// complete. Pass a row to <see cref="DownloadFileAsync(PrideArchiveFile, string, PrideFileChecksum, bool, bool, CancellationToken)"/>
        /// to have the download checked against it.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The list describes the files submitted to the project's root, NOT <see cref="GetProjectFilesAsync"/>'s
        /// manifest, and the two differ in both directions (measured live 2026-09-24). Files PRIDE generated
        /// itself (under <c>generated/</c>) have no row, so a missing row means "cannot be checked", never
        /// "invalid". And the list can carry rows that are not deposits at all: PXD015239's includes a stray
        /// NFS temp file, <c>.nfs80360b91005fae5300002478</c>.
        /// </para>
        /// <para>
        /// PRIDE answers an accession it does not know with 200 and an EMPTY body, exactly as it answers a
        /// project that has no checksums, so a typo would otherwise come back as "nothing to check". On an
        /// empty body this method therefore asks for the project itself (one extra request): if there is no such
        /// project it throws, and only if there is one does it return an empty dictionary. A project not yet
        /// public is unknown to this route too, so it throws as well; PRIDE has no checksums before publication.
        /// </para>
        /// </remarks>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns>
        /// One <see cref="PrideFileChecksum"/> per file, keyed by its bare file name (ordinal, case-sensitive).
        /// Empty only when the project exists and PRIDE holds no checksum for any of its files, in which case
        /// nothing can be verified. Never null.
        /// </returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="MzLibException">
        /// No project has that accession, or PRIDE answered with a list this method cannot read: a header other
        /// than <c>File-Name, File-MD5Checksum, File-Size</c>, an empty line before the last row, a row without exactly those three fields, an MD5 that
        /// is not 32 hexadecimal characters, a size that is not a whole number of bytes, or a file name listed
        /// twice. These are a broken contract, not an outage, so they are not retried.
        /// </exception>
        /// <exception cref="HttpRequestException">The API returned a non-success status code or did not answer. A transient failure is retried, as <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/> describes, before this is thrown.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<IReadOnlyDictionary<string, PrideFileChecksum>> GetFileChecksumsAsync(string accession,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A PRIDE project accession is required.", nameof(accession));

            cancellationToken.ThrowIfCancellationRequested();

            // No Accept header is set, and none may be: this endpoint serves text/plain only and answers
            // "Accept: application/json" with 406.
            string requestUri = $"files/checksum/{Uri.EscapeDataString(accession)}";
            string body;
            using (HttpResponseMessage response = await GetSuccessAsync(requestUri,
                       r => $"PRIDE Archive request failed with status {(int)r.StatusCode} {r.ReasonPhrase} for '{requestUri}'.",
                       cancellationToken).ConfigureAwait(false))
                body = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

            if (string.IsNullOrWhiteSpace(body))
            {
                // An unknown accession and a project without checksums look identical here. GetProjectAsync tells
                // them apart and throws MzLibException for the first, so a typo fails instead of reading as "no rows".
                await GetProjectAsync(accession, cancellationToken).ConfigureAwait(false);
                return new Dictionary<string, PrideFileChecksum>(StringComparer.Ordinal);
            }

            return ParseChecksums(body, accession);
        }

        /// <summary>The header line PRIDE's checksum list starts with, verified live 2026-09-24.</summary>
        internal const string ChecksumHeader = "File-Name\tFile-MD5Checksum\tFile-Size";

        /// <summary>
        /// Reads PRIDE's checksum list strictly: anything other than the known header followed by well-formed
        /// three-field rows is an <see cref="MzLibException"/>, never a silently skipped line, because a row
        /// skipped here would later read as "this file has no checksum" rather than as a broken list.
        /// </summary>
        internal static Dictionary<string, PrideFileChecksum> ParseChecksums(string body, string accession)
        {
            string[] lines = body.Split('\n');
            string header = lines[0].TrimEnd('\r');
            if (header != ChecksumHeader)
                throw new MzLibException(
                    $"PRIDE Archive's checksum list for '{accession}' starts with '{header}', not the expected '{ChecksumHeader.Replace('\t', ',')}'.");

            var checksums = new Dictionary<string, PrideFileChecksum>(StringComparer.Ordinal);
            for (int i = 1; i < lines.Length; i++)
            {
                string line = lines[i].TrimEnd('\r');
                if (line.Length == 0)
                {
                    // Only the newline that ends the last row may leave an empty line; one in the middle means
                    // the list is not what PRIDE has always served.
                    if (lines.Skip(i + 1).All(rest => rest.TrimEnd('\r').Length == 0))
                        break;
                    throw new MzLibException($"PRIDE Archive's checksum list for '{accession}' has an empty line at line {i + 1}.");
                }

                string[] fields = line.Split('\t');
                if (fields.Length != 3 || fields[0].Length == 0)
                    throw new MzLibException(
                        $"PRIDE Archive's checksum list for '{accession}' has a malformed row at line {i + 1}: expected a file name, an MD5 and a size, separated by tabs.");
                if (fields[1].Length != 32 || !fields[1].All(Uri.IsHexDigit))
                    throw new MzLibException(
                        $"PRIDE Archive's checksum list for '{accession}' gives '{fields[0]}' the MD5 '{fields[1]}', which is not 32 hexadecimal characters.");
                if (!long.TryParse(fields[2], NumberStyles.None, CultureInfo.InvariantCulture, out long size))
                    throw new MzLibException(
                        $"PRIDE Archive's checksum list for '{accession}' gives '{fields[0]}' the size '{fields[2]}', which is not a whole number of bytes.");
                if (!checksums.TryAdd(fields[0], new PrideFileChecksum(fields[0], fields[1], size)))
                    throw new MzLibException(
                        $"PRIDE Archive's checksum list for '{accession}' lists '{fields[0]}' twice.");
            }

            return checksums;
        }

        /// <summary>
        /// Finds PRIDE Archive projects matching a free-text keyword (v3 <c>search/projects</c>) — the
        /// discovery entry point for a caller who has a subject rather than an accession.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Hits come back as <see cref="PrideProjectSearchResult"/>, NOT <see cref="PrideProject"/>.
        /// PRIDE serves search from a separate flattened projection in which controlled-vocabulary
        /// terms are reduced to display strings and contacts to names — see the remarks on
        /// <see cref="PrideProjectSearchResult"/>. Follow a hit's
        /// <see cref="PrideProjectSearchResult.Accession"/> to <see cref="GetProjectAsync"/> when the
        /// full metadata object is wanted.
        /// </para>
        /// <para>
        /// Only <paramref name="keyword"/> is exposed. PRIDE also accepts <c>filter</c>,
        /// <c>sortFields</c> and <c>sortDirection</c>, but validates NONE of them: a misspelled field
        /// or an invalid direction returns 200 with unfiltered, unsorted results (verified live
        /// 2026-07-23), so a caller typo would silently produce wrong data instead of an error. Those
        /// parameters are deferred until they can be validated in C# — an enum for the direction, a
        /// restricted set for the sort fields — so a mistake fails here rather than at PRIDE.
        /// </para>
        /// <para>
        /// A keyword is required. PRIDE treats an absent one as "browse the whole archive" (40 000+
        /// projects at the time of writing), which is a different capability with a different cost,
        /// not a degenerate search — so asking for it has to be deliberate rather than the result of
        /// passing through an empty string.
        /// </para>
        /// <para>
        /// EVERY matching project is returned, which for a search means the cost is set by the
        /// KEYWORD rather than by anything the caller can cap. That differs in kind from
        /// <see cref="GetProjectFilesAsync"/>, whose result is bounded by one project. Search hits are
        /// also fat — roughly 10-14 KB each, since each carries the full protocols, every file name and
        /// the match highlights — so a request count badly understates the load. Measured live
        /// 2026-08-21: "liver" was 2 197 projects over 22 requests and about 30 MB; "proteomics" was
        /// 37 356 projects over 374 requests and roughly 355 MB, most of the archive. PRIDE offers no
        /// compression even when asked, so none of that can be traded away. Search narrowly;
        /// <paramref name="pageSize"/> changes only how many requests it takes, never how much comes
        /// back. Narrowing beyond a keyword needs <c>filter</c>, which is deferred for the reason above.
        /// </para>
        /// <para>
        /// The keyword is free text with AND-of-prefix-token semantics, and PRIDE supports no query
        /// operators: quotes are discarded (there is no phrase search), <c>*</c> is dropped, and
        /// <c>AND</c>/<c>OR</c> match as ordinary literal terms — "liver OR kidney" returned 2 hits
        /// where "liver" alone returned 2 197. Diacritics are not folded either: "Nájera" found nothing
        /// while "Najera" found a project. None of this can be escaped around, so a caller who passes
        /// query syntax silently gets a different result set (all verified live 2026-08-21).
        /// </para>
        /// </remarks>
        /// <param name="keyword">
        /// The free-text query. Escaped before it is sent, so it may contain any character —
        /// <c>&amp;</c> and <c>=</c> included, which would otherwise split it into further query
        /// parameters.
        /// </param>
        /// <param name="pageSize">
        /// Hits requested per page (default 100). PRIDE caps this server-side at 100 and then pages by
        /// the capped size, exactly as it does for the file manifest; see
        /// <see cref="GetProjectFilesAsync"/> for what that means for termination. Every page is
        /// fetched regardless, so this is never a limit on how many hits are returned.
        /// <para>
        /// Above the cap it is not a throughput knob either — it is a no-op. A request for 500 comes
        /// back byte-identical to a request for 100 (verified live 2026-08-21), so the fetch costs the
        /// same requests and buys nothing. Below the cap it only costs MORE requests, and a small value
        /// also guarantees the tail probe fires, since every page is then trivially "full". The default
        /// is the value to leave it at.
        /// </para>
        /// </param>
        /// <param name="cancellationToken">Cancels the (possibly multi-page) search.</param>
        /// <returns>
        /// Every matching project, across all pages, with no accession repeated. Empty when nothing
        /// matches — PRIDE reports no hits as an empty result rather than an error, so there is no
        /// <c>Try</c> variant of this method as there is for <see cref="GetProjectAsync"/>. Never null.
        /// <para>
        /// One caveat, stated because it cannot be fixed here: PRIDE pages a LIVE index and offers no
        /// stable cursor, so a result set that changes DURING a multi-page fetch shifts its own paging.
        /// A project published mid-fetch is served on two pages — that is deduplicated, so it comes
        /// back once — but a project REMOVED mid-fetch slides the window the other way and can fall
        /// between two pages, and then no page carries it. A search whose results fit on one page
        /// cannot be affected. Verified live 2026-08-21.
        /// </para>
        /// </returns>
        /// <exception cref="ArgumentException">
        /// The keyword is null, empty, whitespace, or longer than <see cref="MaxKeywordLength"/>.
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException">The page size is not positive.</exception>
        /// <exception cref="HttpRequestException">The API returned a non-success status code or did not answer. A transient failure of one page is retried, that page alone, as <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/> describes before this is thrown.</exception>
        /// <exception cref="MzLibException">
        /// PRIDE answered successfully but served a page identical to its predecessor while
        /// <c>total_records</c> reported more remained — a broken contract rather than an outage.
        /// </exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<List<PrideProjectSearchResult>> SearchProjectsAsync(string keyword, int pageSize = 100,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(keyword))
                throw new ArgumentException("A search keyword is required.", nameof(keyword));
            if (keyword.Length > MaxKeywordLength)
                throw new ArgumentException(
                    $"A search keyword may be at most {MaxKeywordLength} characters, but this one is {keyword.Length}. " +
                    "PRIDE answers a very long keyword with HTTP 500, which is indistinguishable from the service " +
                    "being down.", nameof(keyword));
            if (pageSize <= 0)
                throw new ArgumentOutOfRangeException(nameof(pageSize), "Page size must be positive.");

            List<PrideProjectSearchResult> hits = await GetAllPagesAsync<PrideProjectSearchResult>(
                page => $"search/projects?keyword={Uri.EscapeDataString(keyword)}&pageSize={pageSize}&page={page}",
                pageSize,
                $"keyword '{keyword}'",
                cancellationToken).ConfigureAwait(false);

            // PRIDE pages a LIVE index with no stable cursor, so the result set can change underneath a
            // multi-page fetch. Verified live 2026-08-21: a project left the "liver" set mid-session and
            // every record after it shifted up one position, which moves a record from page 1 into
            // page 0's range after page 0 was already read. Publication does the same in reverse, and
            // then a record is served on two consecutive pages.
            //
            // The shared pager cannot fix this, and deliberately does not try: it refuses to treat any
            // FIELD as a record's identity because the file manifest it was written for has no unique
            // key, so dropping "repeats" there would drop real files. That reasoning does not carry
            // over. A search hit HAS a real identity — its accession is the very thing a caller would
            // use to fetch the project — so two hits sharing one are the same project, not two projects
            // that happen to look alike. Deduping here is therefore sound where it would not be there.
            //
            // Only the duplicate half is fixable. A record can equally be skipped, when the window
            // slides the other way, and nothing short of a cursor the API does not offer would catch
            // that. It is stated in the returns docs rather than papered over.
            var seenAccessions = new HashSet<string>(StringComparer.Ordinal);
            return hits.Where(hit => seenAccessions.Add(hit.Accession)).ToList();
        }

        /// <summary>
        /// Downloads a single PRIDE file's bytes to <paramref name="destinationDirectory"/>, saved under the
        /// file's own <see cref="PrideArchiveFile.FileName"/>. The download runs over HTTPS through this
        /// client's reused <see cref="HttpClient"/>: PRIDE exposes files as FTP/Aspera locations, but its FTP
        /// host also serves the identical path over HTTPS, so an FTP location is scheme-upgraded to HTTPS (see
        /// <see cref="PrideArchiveExtensions.GetHttpsDownloadUrl"/>). The bytes are streamed to a sibling
        /// ".partial" file and moved into place only on success, so an interrupted transfer never leaves a
        /// truncated file at the destination path.
        /// </summary>
        /// <remarks>
        /// A transient failure is retried up to three times, after 5 s, 20 s and 60 s (a 429 or 503 that sends
        /// <c>Retry-After</c> sets its own wait, capped at 60 s). Transient means a timeout, a dropped or stalled
        /// body, 408, 429, any 5xx, and 403 from PRIDE's FTP host, which is how EBI rate-limits it.
        /// <para>
        /// A retry resumes rather than restarting: it asks only for the bytes the ".partial" lacks, under
        /// <c>If-Range</c>, so a file that changed on the server is fetched whole instead of spliced. When the
        /// retries run out the ".partial" is KEPT, beside a ".partial.validator" file recording the validator it
        /// was served under, and the next call for the same file resumes it. That is what keeps a re-run after a
        /// failed multi-gigabyte transfer cheap. A ".partial" with no validator (a server that sends neither
        /// <c>ETag</c> nor <c>Last-Modified</c>) cannot be resumed safely, so it is not kept. After any failure that
        /// is not transient, both files are deleted.
        /// </para>
        /// </remarks>
        /// <param name="file">The file to download. Must not be null and must have a file name.</param>
        /// <param name="destinationDirectory">The directory to write into; created if it does not exist.</param>
        /// <param name="overwrite">
        /// When true (the default) an existing destination file is replaced. When false, an existing
        /// destination file is left untouched and no request is made (a cheap resume for large projects).
        /// </param>
        /// <param name="cancellationToken">Cancels the download.</param>
        /// <returns>The full path of the written (or already-present) file.</returns>
        /// <exception cref="ArgumentNullException">The file is null.</exception>
        /// <exception cref="ArgumentException">The destination directory is blank, the file has no name, or the file name is not a bare file name (contains a path separator, a "..", or a root).</exception>
        /// <exception cref="NotSupportedException">The file exposes no HTTPS-reachable location (e.g. Aspera-only).</exception>
        /// <exception cref="HttpRequestException">The download returned a non-success status code (carried in <see cref="HttpRequestException.StatusCode"/>), PRIDE did not answer within the client's timeout, the connection dropped while the body was being read, the body delivered nothing for <see cref="BodyStallTimeout"/>, or the file did not come to the length the server announced. A transient one is thrown only after its retries are spent, as the last attempt's exception. The message names the file and host, never the URL.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled. A stall is NOT reported this way — see <see cref="BodyStallTimeout"/>.</exception>
        public Task<string> DownloadFileAsync(PrideArchiveFile file, string destinationDirectory,
            bool overwrite = true, CancellationToken cancellationToken = default) =>
            DownloadFileCoreAsync(file, destinationDirectory, expected: null, overwrite, verifyMd5: false, cancellationToken);

        /// <summary>
        /// Downloads a single PRIDE file exactly as <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/>
        /// does, and checks it against the size, and optionally the MD5, that PRIDE recorded for it
        /// (<see cref="GetFileChecksumsAsync"/>). A file that fails the check never reaches the destination path.
        /// </summary>
        /// <remarks>
        /// The size is always checked. The bytes are checked too, as they are written, so neither check reads the
        /// file a second time: with <paramref name="verifyMd5"/> true against PRIDE's MD5, and otherwise by a scan
        /// for the error response EBI's storage has been seen writing into files it serves at their full size
        /// (<c>&lt;Error&gt;&lt;Code&gt;</c>, 2026-10-08), which a size check cannot see. A resumed transfer is read
        /// once more at the end, because its start arrived in an earlier attempt.
        /// <para>
        /// With <paramref name="overwrite"/> false, an existing destination file is kept only if it passes the size
        /// check, and the MD5 when <paramref name="verifyMd5"/> is true (reading it again, which takes seconds for a
        /// multi-gigabyte raw file). It is not scanned. One that fails (a truncated copy from an older tool, say) is
        /// downloaded again and replaced, instead of being returned as if it were complete.
        /// </para>
        /// <para>
        /// A download that fails the check is deleted, with its ".partial.validator", and downloaded again from its
        /// first byte, never resumed, within the same retry budget as a transient failure. EBI's damage is
        /// transient: a fresh download of the same file usually comes out right. Only when every attempt fails does
        /// it throw <see cref="MzLibException"/>, leaving no ".partial" behind.
        /// </para>
        /// <para>
        /// Cancelling while the downloaded file is being checked keeps its ".partial" and ".partial.validator".
        /// The next call for the same file with a checksum whose size the ".partial" already has checks it again
        /// instead of downloading it again, so cancelling the MD5 of a multi-gigabyte file does not throw the
        /// transfer away.
        /// </para>
        /// </remarks>
        /// <param name="file">The file to download. Must not be null and must have a file name.</param>
        /// <param name="destinationDirectory">The directory to write into; created if it does not exist.</param>
        /// <param name="expected">PRIDE's checksum row for this file. Its <see cref="PrideFileChecksum.FileName"/> must equal <paramref name="file"/>'s.</param>
        /// <param name="overwrite">
        /// When true (the default) an existing destination file is replaced. When false, an existing destination
        /// file that passes the check is left untouched and no request is made.
        /// </param>
        /// <param name="verifyMd5">When true, the MD5 is checked as well as the size. False by default.</param>
        /// <param name="cancellationToken">Cancels the download, or the hashing of an existing file.</param>
        /// <returns>The full path of the written (or already-present and checked) file.</returns>
        /// <exception cref="ArgumentNullException">The file or <paramref name="expected"/> is null.</exception>
        /// <exception cref="ArgumentException">As for the overload without a checksum, or <paramref name="expected"/> is for a different file name.</exception>
        /// <exception cref="NotSupportedException">The file exposes no HTTPS-reachable location (e.g. Aspera-only) and must be downloaded.</exception>
        /// <exception cref="HttpRequestException">As for the overload without a checksum.</exception>
        /// <exception cref="MzLibException">On every attempt, the downloaded file's size differed from <paramref name="expected"/>, or its MD5 did when <paramref name="verifyMd5"/> is true, or (without the MD5) it held EBI's error text. The message, the last attempt's, names the file and host, never the URL.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled.</exception>
        public Task<string> DownloadFileAsync(PrideArchiveFile file, string destinationDirectory,
            PrideFileChecksum expected, bool overwrite = true, bool verifyMd5 = false,
            CancellationToken cancellationToken = default)
        {
            if (expected == null)
                throw new ArgumentNullException(nameof(expected));
            return DownloadFileCoreAsync(file, destinationDirectory, expected, overwrite, verifyMd5, cancellationToken);
        }

        private async Task<string> DownloadFileCoreAsync(PrideArchiveFile file, string destinationDirectory,
            PrideFileChecksum expected, bool overwrite, bool verifyMd5, CancellationToken cancellationToken)
        {
            if (file == null)
                throw new ArgumentNullException(nameof(file));
            if (string.IsNullOrWhiteSpace(destinationDirectory))
                throw new ArgumentException("A destination directory is required.", nameof(destinationDirectory));
            if (string.IsNullOrWhiteSpace(file.FileName))
                throw new ArgumentException("The PRIDE file has no file name to save under.", nameof(file));

            // The file name comes verbatim from the PRIDE response; treat it as untrusted. Only a bare
            // leaf name is allowed, so a value carrying a directory separator, a ".." segment, or a rooted
            // path cannot escape destinationDirectory when combined below (Path.Combine does not sanitize).
            string safeFileName = Path.GetFileName(file.FileName);
            if (safeFileName != file.FileName || safeFileName == "." || safeFileName == "..")
                throw new ArgumentException(
                    $"The PRIDE file name '{file.FileName}' is not a bare file name; refusing to write outside the destination directory.",
                    nameof(file));

            if (expected != null && !string.Equals(expected.FileName, safeFileName, StringComparison.Ordinal))
                throw new ArgumentException(
                    $"The checksum given is for '{expected.FileName}', not for '{safeFileName}'.", nameof(expected));

            string destinationPath = Path.Combine(destinationDirectory, safeFileName);

            // Cheap resume: an already-present destination is left untouched. This runs before URL
            // resolution and directory creation so skipping a downloaded file never fails on a file
            // that has no HTTPS location (e.g. Aspera-only) or does needless filesystem work. With a
            // checksum, only a file that passes it counts as already present.
            if (!overwrite && File.Exists(destinationPath)
                && (expected == null
                    || await ChecksumMismatchAsync(destinationPath, expected, verifyMd5, cancellationToken).ConfigureAwait(false) == null))
                return destinationPath;

            string url = file.GetHttpsDownloadUrl(); // throws NotSupportedException if unreachable over HTTPS
            string described = DescribeDownload(url, safeFileName);

            Directory.CreateDirectory(destinationDirectory);

            string partialPath = destinationPath + ".partial";
            string validatorPath = partialPath + ValidatorSuffix;
            string host = HostOf(url);

            // A .partial left by an earlier call is resumable only with the validator it was started under;
            // without one, nothing proves the bytes on disk are the start of the file the server holds now.
            var transfer = new ResumableTransfer
            {
                Validator = ReadValidator(partialPath, validatorPath),
                Digest = expected != null ? new DownloadDigest(verifyMd5) : null
            };
            bool verifying = false;

            try
            {
                // A .partial that already has the listed size was downloaded whole under its validator by a call
                // cancelled while checking it (see the catch below): check it again rather than fetch it again.
                // If it fails, it is discarded and the file is downloaded as if it had never been there.
                bool complete = false;
                if (expected != null && transfer.Validator != null && new FileInfo(partialPath).Length == expected.SizeBytes)
                {
                    verifying = true;
                    transfer.Digest.Invalidate();
                    complete = await ContentMismatchAsync(partialPath, expected, verifyMd5, transfer.Digest, cancellationToken).ConfigureAwait(false) == null;
                    if (!complete)
                        Discard(partialPath, validatorPath, transfer);
                    verifying = false;
                }

                if (!complete)
                {
                    await WithRetryAsync(async () =>
                    {
                        verifying = false;
                        await DownloadOnceAsync(url, described, partialPath, validatorPath, transfer, cancellationToken).ConfigureAwait(false);
                        if (expected == null)
                            return true;

                        // Checked before the move, so a file that fails never reaches the destination path. The
                        // check is made inside the attempt so that a failure is retried: EBI has served files at
                        // their full size with its own storage error text in place of some of their bytes, and a
                        // fresh download of the same file usually comes out right (measured 2026-10-08).
                        verifying = true;
                        cancellationToken.ThrowIfCancellationRequested();
                        string mismatch = await ContentMismatchAsync(partialPath, expected, verifyMd5, transfer.Digest, cancellationToken).ConfigureAwait(false);
                        if (mismatch != null)
                        {
                            // Never resumed: the bytes on disk are what is wrong.
                            Discard(partialPath, validatorPath, transfer);
                            throw new ContentMismatchException($"The PRIDE download of {described} {mismatch}.");
                        }
                        return true;
                    }, host, cancellationToken).ConfigureAwait(false);
                }

                File.Move(partialPath, destinationPath, overwrite: true);
                DeleteQuietly(validatorPath);
            }
            catch (ContentMismatchException e)
            {
                // Every attempt came back wrong, so this is no longer an outage worth waiting out: PRIDE's file and
                // PRIDE's own list disagree. Nothing is kept; the partial went with the attempt that failed.
                DeleteQuietly(partialPath);
                DeleteQuietly(validatorPath);
                throw new MzLibException(e.Message, e);
            }
            catch (HttpRequestException e) when (IsTransient(e, host) && transfer.Validator != null && File.Exists(partialPath))
            {
                // Kept on purpose: EBI failed, not the file, and the validator beside it lets the next call
                // ask only for the bytes that are missing instead of paying for the whole file again.
                throw;
            }
            catch (OperationCanceledException) when (verifying && transfer.Validator != null && File.Exists(partialPath))
            {
                // Kept on purpose: the caller stopped the check, not the transfer, and the whole file is on disk.
                // Deleting it would throw away a complete multi-gigabyte download.
                throw;
            }
            catch
            {
                // Anything else -- a refusal, a broken contract, the caller's disk, the caller's cancellation --
                // gives no reason to trust the partial bytes, so neither they nor their validator survive.
                DeleteQuietly(partialPath);
                DeleteQuietly(validatorPath);
                throw;
            }
            finally
            {
                transfer.Digest?.Dispose();
            }

            return destinationPath;
        }

        /// <summary>
        /// The suffix of the sidecar file, beside a <c>.partial</c>, that holds the validator (strong
        /// <c>ETag</c>, else <c>Last-Modified</c>) the partial bytes were downloaded under.
        /// </summary>
        internal const string ValidatorSuffix = ".validator";

        /// <summary>
        /// Checks a file already at the destination against PRIDE's checksum row: its length always, its MD5 only
        /// when <paramref name="verifyMd5"/> is set. Returns null when it passes, else the end of a sentence saying
        /// what differs (it names no path or URL, so it can go into an exception message as it is).
        /// </summary>
        /// <remarks>
        /// Unlike a fresh download, which is scanned for EBI's error text as it arrives, a file already on disk is
        /// not read at all unless its MD5 is asked for, so skipping a project's thousands of finished files stays cheap.
        /// </remarks>
        private static async Task<string> ChecksumMismatchAsync(string path, PrideFileChecksum expected, bool verifyMd5,
            CancellationToken cancellationToken)
        {
            long length = new FileInfo(path).Length;
            if (length != expected.SizeBytes)
                return SizeMismatch(length, expected);
            if (!verifyMd5)
                return null;

            var digest = new DownloadDigest(verifyMd5: true);
            digest.Invalidate();
            return await ContentMismatchAsync(path, expected, verifyMd5: true, digest, cancellationToken).ConfigureAwait(false);
        }

        /// <summary>
        /// Checks a downloaded file against PRIDE's checksum row: its length always, then its MD5 when
        /// <paramref name="verifyMd5"/> is set, else a scan for <see cref="DownloadDigest.ServerErrorMarker"/>. Uses
        /// what <paramref name="digest"/> saw as the bytes arrived, and reads the file again only when the digest did
        /// not see all of it (a resumed transfer, or a partial left by an earlier call). Returns null when it passes,
        /// else the end of a sentence saying what differs, naming no path or URL.
        /// </summary>
        private static async Task<string> ContentMismatchAsync(string path, PrideFileChecksum expected, bool verifyMd5,
            DownloadDigest digest, CancellationToken cancellationToken)
        {
            long length = new FileInfo(path).Length;
            if (length != expected.SizeBytes)
                return SizeMismatch(length, expected);

            if (!digest.SawWholeFile)
            {
                digest.Reset();
                byte[] buffer = new byte[1 << 16];
                using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.Read, buffer.Length, useAsync: true);
                int read;
                while ((read = await stream.ReadAsync(buffer.AsMemory(), cancellationToken).ConfigureAwait(false)) > 0)
                    digest.Append(buffer.AsSpan(0, read));
            }

            if (verifyMd5)
            {
                // A matching MD5 proves the file is PRIDE's, so the error-text scan has nothing to add.
                string md5 = digest.Md5Hex();
                return string.Equals(md5, expected.Md5, StringComparison.OrdinalIgnoreCase)
                    ? null
                    : $"has MD5 {md5} where PRIDE's checksum list records {expected.Md5}";
            }

            return digest.ServerErrorAt >= 0
                ? $"holds a server error response (\"<Error><Code>\") at byte {digest.ServerErrorAt}, written into the file in place of its own bytes"
                : null;
        }

        private static string SizeMismatch(long length, PrideFileChecksum expected) =>
            $"came to {length} bytes where PRIDE's checksum list records {expected.SizeBytes}";

        /// <summary>
        /// Deletes a partial whose bytes are wrong, with its validator, so the next attempt starts from zero instead
        /// of resuming them.
        /// </summary>
        private static void Discard(string partialPath, string validatorPath, ResumableTransfer transfer)
        {
            DeleteQuietly(partialPath);
            DeleteQuietly(validatorPath);
            transfer.Validator = null;
        }

        /// <summary>What one download carries from a failed attempt into the next.</summary>
        private sealed class ResumableTransfer
        {
            /// <summary>The validator the partial bytes were served under, or null when they cannot be resumed.</summary>
            public RangeConditionHeaderValue Validator { get; set; }

            /// <summary>The check run over the bytes as they arrive, or null when there is nothing to check them against.</summary>
            public DownloadDigest Digest { get; init; }
        }

        /// <summary>
        /// A downloaded file is not the attempt's failure but its content's: wrong size, wrong MD5, or EBI's error text
        /// inside it. It has no status, so <see cref="WithRetryAsync{T}"/> treats it as transient and downloads the file
        /// again from zero; once the retries are spent it becomes an <see cref="MzLibException"/>.
        /// </summary>
        private sealed class ContentMismatchException(string message) : HttpRequestException(message);

        /// <summary>
        /// The MD5 of a download, or a scan of it for <see cref="ServerErrorMarker"/>, taken in the same pass that
        /// writes it, so checking a multi-gigabyte file costs no second read.
        /// </summary>
        /// <remarks>
        /// EBI's storage has served files at their full, correct size with blocks of its own error response written
        /// over some of their bytes: <c>&lt;Error&gt;&lt;Code&gt;ConnectionClosedException&lt;/Code&gt;&lt;Message&gt;Premature
        /// end of Content-Length delimited message body ...&lt;/Message&gt;...&lt;/Error&gt;</c>, one to four times per file,
        /// over HTTPS and Aspera alike, while every client reported success (PXReprise, 2026-10-07/08; 15 damaged
        /// copies kept). A size check cannot see it. The MD5 can, when the caller asks for it; the scan catches it
        /// without one. A raw file has no reason to hold that text.
        /// </remarks>
        internal sealed class DownloadDigest : IDisposable
        {
            /// <summary>The start of every error response EBI's storage wrote into served files.</summary>
            internal static readonly byte[] ServerErrorMarker = "<Error><Code>"u8.ToArray();

            private readonly bool _verifyMd5;
            private readonly byte[] _tail = new byte[ServerErrorMarker.Length - 1];
            private IncrementalHash _md5;
            private int _tailLength;
            private long _length;

            /// <param name="verifyMd5">Take the MD5; otherwise scan for <see cref="ServerErrorMarker"/>.</param>
            internal DownloadDigest(bool verifyMd5)
            {
                _verifyMd5 = verifyMd5;
                Reset();
            }

            /// <summary>Whether every byte of the file went through <see cref="Append"/> since the last <see cref="Reset"/>.</summary>
            internal bool SawWholeFile { get; private set; }

            /// <summary>The offset of the first <see cref="ServerErrorMarker"/> seen, or -1.</summary>
            internal long ServerErrorAt { get; private set; }

            /// <summary>Starts over, for a file that is about to arrive from its first byte.</summary>
            internal void Reset()
            {
                _md5?.Dispose();
                _md5 = _verifyMd5 ? IncrementalHash.CreateHash(HashAlgorithmName.MD5) : null;
                _tailLength = 0;
                _length = 0;
                ServerErrorAt = -1;
                SawWholeFile = true;
            }

            /// <summary>Marks the file as not seen whole (its start is already on disk), so it must be read again to be checked.</summary>
            internal void Invalidate() => SawWholeFile = false;

            /// <summary>Takes in the next bytes of the file, in order.</summary>
            internal void Append(ReadOnlySpan<byte> data)
            {
                if (_md5 != null)
                    _md5.AppendData(data);
                else if (ServerErrorAt < 0)
                    Scan(data);
                _length += data.Length;
            }

            public void Dispose() => _md5?.Dispose();

            /// <summary>The MD5 of everything appended, as lowercase hex.</summary>
            internal string Md5Hex() => Convert.ToHexString(_md5.GetHashAndReset()).ToLowerInvariant();

            /// <summary>
            /// Looks for the marker in <paramref name="data"/>, and across the boundary with the previous read: the
            /// last <c>marker length - 1</c> bytes are carried over, as <c>MzmlMethods.GetSHA1Hash</c> does for its tag.
            /// </summary>
            private void Scan(ReadOnlySpan<byte> data)
            {
                int carry = _tail.Length;
                Span<byte> seam = stackalloc byte[2 * carry];
                _tail.AsSpan(0, _tailLength).CopyTo(seam);
                int head = Math.Min(carry, data.Length);
                data[..head].CopyTo(seam[_tailLength..]);
                int at = seam[..(_tailLength + head)].IndexOf(ServerErrorMarker);
                if (at >= 0)
                {
                    ServerErrorAt = _length - _tailLength + at;
                    return;
                }

                at = data.IndexOf(ServerErrorMarker);
                if (at >= 0)
                {
                    ServerErrorAt = _length + at;
                    return;
                }

                if (data.Length >= carry)
                {
                    data[^carry..].CopyTo(_tail);
                    _tailLength = carry;
                }
                else
                {
                    // Keep the last `carry` bytes of what was carried plus this short read.
                    int keep = Math.Min(carry, _tailLength + data.Length);
                    seam[(_tailLength + data.Length - keep)..(_tailLength + data.Length)].CopyTo(_tail);
                    _tailLength = keep;
                }
            }
        }

        /// <summary>
        /// One download attempt into <paramref name="partialPath"/>: a ranged request when the partial can be
        /// resumed, else a full one. Returns once the partial holds the complete file.
        /// </summary>
        /// <remarks>
        /// Bytes are appended only to a <c>206</c> whose <c>Content-Range</c> starts exactly where the partial
        /// ends, sent under <c>If-Range</c> so the server answers with the whole file instead if it has changed.
        /// Any other answer to a ranged request -- a <c>200</c>, a <c>206</c> starting elsewhere, a <c>416</c> --
        /// starts the file again from zero. The guarantee is that no truncated or spliced file is ever moved
        /// into place, not that no <c>.partial</c> is ever left behind.
        /// </remarks>
        private async Task<bool> DownloadOnceAsync(string url, string described, string partialPath,
            string validatorPath, ResumableTransfer transfer, CancellationToken cancellationToken)
        {
            long resumeFrom = transfer.Validator != null && File.Exists(partialPath) ? new FileInfo(partialPath).Length : 0;

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
                        $"PRIDE download failed with status {(int)response.StatusCode} {response.ReasonPhrase} for {described}.");

                long? expectedLength;
                if (resumed)
                {
                    expectedLength = response.Content.Headers.ContentRange.Length;
                    transfer.Digest?.Invalidate(); // the start of the file is on disk, not in this stream
                }
                else
                {
                    // A fresh start: whatever is on disk is discarded, and the validator is the new response's.
                    expectedLength = response.Content.Headers.ContentLength;
                    transfer.Validator = ValidatorOf(response);
                    WriteValidator(validatorPath, transfer.Validator);
                    transfer.Digest?.Reset();
                }

                using (var fileStream = new FileStream(partialPath, resumed ? FileMode.Append : FileMode.Create,
                           FileAccess.Write, FileShare.None))
                using (Stream httpStream = await response.Content.ReadAsStreamAsync(cancellationToken).ConfigureAwait(false))
                {
                    // Not Stream.CopyToAsync: it would inherit the very absence of a read deadline that
                    // BodyStallTimeout exists to supply, which is how the body escaped every timeout here.
                    await CopyUntilStalledAsync(httpStream, fileStream, described, transfer.Digest, cancellationToken).ConfigureAwait(false);
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
                        $"The PRIDE download of {described} came to {actualLength} bytes where the server announced {expectedLength.Value}.");
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
        /// and only if it still matches <paramref name="validator"/>.
        /// </summary>
        private Task<HttpResponseMessage> SendDownloadRequestAsync(string url, string described, long resumeFrom,
            RangeConditionHeaderValue validator, CancellationToken cancellationToken)
        {
            var request = new HttpRequestMessage(HttpMethod.Get, url);
            if (resumeFrom > 0 && validator != null)
            {
                request.Headers.Range = new RangeHeaderValue(resumeFrom, null);
                request.Headers.IfRange = validator;
            }

            return SendAsync(request, HttpCompletionOption.ResponseHeadersRead, cancellationToken, described);
        }

        /// <summary>
        /// The validator to resume a response's body under: its strong <c>ETag</c>, else its <c>Last-Modified</c>,
        /// else null. A weak ETag is not allowed in <c>If-Range</c>.
        /// </summary>
        /// <remarks>
        /// With no validator the bytes cannot be resumed, not even within this call: a <c>Range</c> request with
        /// no <c>If-Range</c> could splice two versions of a file together. PRIDE's public FTP host sends both
        /// validators. Its reviewer-token route sends neither (measured 2026-09-29), so a private download
        /// restarts from zero after a failure.
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
        /// validator, or a validator with no partial -- is deleted, because it can never be resumed.
        /// </summary>
        private static RangeConditionHeaderValue ReadValidator(string partialPath, string validatorPath)
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
        private static void DeleteQuietly(string path)
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
        /// Sends a GET, reporting <see cref="HttpClient.Timeout"/> expiring as the transport failure it is.
        /// </summary>
        /// <remarks>
        /// When the client's own timeout fires, <see cref="HttpClient"/> throws a
        /// <see cref="TaskCanceledException"/> -- the same type a caller's cancellation produces -- so an EBI
        /// that never answers escaped the documented contract (and <c>ExternalServiceTestHelper</c>, which
        /// skips only on transport failures). HttpClient marks its own timeout with an inner
        /// <see cref="TimeoutException"/>, so only that is converted; a caller's cancellation, through the token or
        /// through <see cref="HttpClient.CancelPendingRequests"/> on an injected client, is not.
        /// </remarks>
        /// <param name="requestUri">The URI to fetch.</param>
        /// <param name="completionOption">When the returned task completes.</param>
        /// <param name="cancellationToken">The caller's token; its cancellation is never converted.</param>
        /// <param name="described">How the request is named in the message; defaults to the URI itself.</param>
        private Task<HttpResponseMessage> GetAsync(string requestUri, HttpCompletionOption completionOption,
            CancellationToken cancellationToken, string described = null) =>
            SendAsync(new HttpRequestMessage(HttpMethod.Get, requestUri), completionOption, cancellationToken,
                described ?? $"'{requestUri}'");

        /// <summary>
        /// <see cref="GetAsync"/> for a request that carries headers of its own (a resumed download's
        /// <c>Range</c>). Takes ownership of <paramref name="request"/>.
        /// </summary>
        private async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, HttpCompletionOption completionOption,
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
                        $"PRIDE did not respond within {_httpClient.Timeout} for {described}.", e);
                }
            }
        }

        /// <summary>
        /// Fetches <paramref name="requestUri"/> in full, retrying transient failures, and returns the response
        /// only once it has a success status. Any other status throws, worded by <paramref name="describeFailure"/>.
        /// </summary>
        private Task<HttpResponseMessage> GetSuccessAsync(string requestUri,
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
        private static HttpRequestException StatusFailure(HttpResponseMessage response, string message)
        {
            var exception = new HttpRequestException(message, null, response.StatusCode);
            TimeSpan? retryAfter = response.Headers.RetryAfter?.Delta
                ?? (response.Headers.RetryAfter?.Date is DateTimeOffset date ? date - DateTimeOffset.UtcNow : null);
            if (retryAfter.HasValue)
                exception.Data[RetryAfterKey] = retryAfter.Value;
            return exception;
        }

        private const string RetryAfterKey = "PrideArchiveClient.RetryAfter";

        /// <summary>
        /// Runs <paramref name="attempt"/>, retrying it up to <see cref="MaxRetries"/> times while it fails
        /// transiently, and lets the last failure propagate unchanged: same type, same
        /// <see cref="HttpRequestException.StatusCode"/>, same message.
        /// </summary>
        /// <remarks>
        /// PR-A made every transport failure an <see cref="HttpRequestException"/>, so the decision is made on
        /// the exception alone. One with no status is a transport failure (a timeout, a dropped or stalled
        /// body) and is always transient. One with a status is transient for 408, 429 and every 5xx -- and for
        /// 403 on the FTP host only, because EBI answers rate limiting there with 403 rather than 429. On the
        /// REST host a 403 is a real refusal, and retrying it would only delay the answer by 85 s.
        /// <para>
        /// Never retried: a <see cref="MzLibException"/> (PRIDE answered, and broke its contract, which a
        /// second ask will not mend), a write-side <see cref="IOException"/> (the caller's disk) and the
        /// caller's cancellation, which also ends a backoff in progress.
        /// </para>
        /// </remarks>
        /// <param name="attempt">One complete attempt, from sending the request to checking the status.</param>
        /// <param name="host">The host the request goes to; decides whether a 403 is transient.</param>
        /// <param name="cancellationToken">Cancels the attempt and any backoff.</param>
        private async Task<T> WithRetryAsync<T>(Func<Task<T>> attempt, string host, CancellationToken cancellationToken)
        {
            for (int retry = 0; ; retry++)
            {
                try
                {
                    return await attempt().ConfigureAwait(false);
                }
                catch (HttpRequestException e) when (retry < MaxRetries && IsTransient(e, host))
                {
                    await RetryDelay(BackoffFor(retry, e), cancellationToken).ConfigureAwait(false);
                }
            }
        }

        /// <summary>Whether <see cref="WithRetryAsync{T}"/> would retry <paramref name="exception"/> against <paramref name="host"/>.</summary>
        private static bool IsTransient(HttpRequestException exception, string host)
        {
            if (exception.StatusCode is not HttpStatusCode status)
                return true;

            int code = (int)status;
            return status == HttpStatusCode.RequestTimeout
                || status == HttpStatusCode.TooManyRequests
                || code >= 500
                || (status == HttpStatusCode.Forbidden
                    && string.Equals(host, PrideArchiveExtensions.PrideFtpHost, StringComparison.OrdinalIgnoreCase));
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
        private string HostOf(string requestUri) =>
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
        /// </remarks>
        private async Task CopyUntilStalledAsync(Stream source, Stream destination, string described,
            DownloadDigest digest, CancellationToken cancellationToken)
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
                            $"The PRIDE response body for {described} delivered nothing for {BodyStallTimeout}.", e);
                    }
                    catch (IOException e)
                    {
                        // A connection EBI drops mid-body arrives here as an IOException (HttpIOException on
                        // .NET 10, "The response ended prematurely"), not as the HttpRequestException the
                        // contract promises -- so a caller catching "try again later" missed the phase where a
                        // large download actually breaks, and a live test reddened instead of skipping. Only
                        // the READ is wrapped: a failure writing the local file is the caller's disk, not an
                        // outage, and stays an IOException.
                        throw new HttpRequestException(
                            $"The PRIDE response body for {described} ended before the download completed.", e);
                    }
                }

                if (read == 0)
                    return;

                digest?.Append(buffer.AsSpan(0, read));
                await destination.WriteAsync(buffer.AsMemory(0, read), cancellationToken).ConfigureAwait(false);
            }
        }

        /// <summary>
        /// Downloads a project's files to <paramref name="destinationDirectory"/>, optionally filtered. This is
        /// the convenience over <see cref="GetProjectFilesAsync"/> + <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/>: it fetches
        /// the manifest, applies <paramref name="filter"/>, and downloads each selected file in turn.
        /// </summary>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="destinationDirectory">The directory to write into; created if it does not exist.</param>
        /// <param name="filter">
        /// An optional predicate selecting which files to download (e.g. by category or extension — see
        /// <see cref="PrideArchiveExtensions"/>). When null, every file in the manifest is downloaded.
        /// </param>
        /// <param name="overwrite">Passed through to <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/>; default true.</param>
        /// <param name="cancellationToken">Cancels between and during file downloads.</param>
        /// <returns>The full paths of the downloaded files, in manifest order. Empty if none matched.</returns>
        /// <exception cref="ArgumentException">The accession or destination directory is blank.</exception>
        /// <exception cref="HttpRequestException">The manifest request or a download returned a non-success status.</exception>
        public async Task<IReadOnlyList<string>> DownloadProjectFilesAsync(string accession, string destinationDirectory,
            Func<PrideArchiveFile, bool> filter = null, bool overwrite = true, CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(destinationDirectory))
                throw new ArgumentException("A destination directory is required.", nameof(destinationDirectory));

            List<PrideArchiveFile> files = await GetProjectFilesAsync(accession, cancellationToken: cancellationToken).ConfigureAwait(false);
            IEnumerable<PrideArchiveFile> selected = filter == null ? files : files.Where(filter);

            var downloadedPaths = new List<string>();
            foreach (PrideArchiveFile file in selected)
            {
                cancellationToken.ThrowIfCancellationRequested();
                downloadedPaths.Add(await DownloadFileAsync(file, destinationDirectory, overwrite, cancellationToken).ConfigureAwait(false));
            }
            return downloadedPaths;
        }

        /// <summary>
        /// Fetches a single spectrum by its USI (Universal Spectrum Identifier) from PRIDE's PROXI API, returning
        /// the raw PROXI object: the peak arrays plus the controlled-vocabulary
        /// <see cref="PrideProxiSpectrum.Attributes"/> (charge, precursor m/z, ms level, scan number,
        /// instrument, ...). Use <see cref="GetSpectrumAsync"/> instead if you only need the peaks as an
        /// <see cref="MzSpectrum"/> and can discard the attributes.
        /// </summary>
        /// <param name="usi">
        /// The Universal Spectrum Identifier, e.g.
        /// "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2".
        /// </param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns>The spectrum identified by <paramref name="usi"/>. Never null.</returns>
        /// <exception cref="ArgumentException">The USI is null, empty, or whitespace.</exception>
        /// <exception cref="HttpRequestException">
        /// The API returned a non-success status — PROXI answers an unknown or unreadable USI with 404 and a
        /// malformed USI with 400 — or returned an empty result for the USI.
        /// </exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<PrideProxiSpectrum> GetProxiSpectrumAsync(string usi, CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(usi))
                throw new ArgumentException("A USI (Universal Spectrum Identifier) is required.", nameof(usi));

            cancellationToken.ThrowIfCancellationRequested();

            // PROXI is a different path root than the v3 archive BaseAddress, so an archive-relative URI (the
            // GetProjectFilesAsync pattern) would resolve to the wrong path. An absolute PROXI URI overrides the
            // client's BaseAddress. resultType=full asks PROXI for the peak arrays, not just the metadata.
            string requestUri = $"{DefaultProxiBaseAddress}spectra?usi={Uri.EscapeDataString(usi)}&resultType=full";
            using HttpResponseMessage response = await GetSuccessAsync(requestUri,
                r => $"PRIDE PROXI request failed with status {(int)r.StatusCode} {r.ReasonPhrase} for '{requestUri}'.",
                cancellationToken).ConfigureAwait(false);

            string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

            // PROXI wraps its result in a JSON array (one entry per matched spectrum). A USI identifies exactly
            // one spectrum, so take the first. An empty array on a 200 is a contract oddity — unknown USIs 404
            // rather than returning [] — so treat "no spectrum" as an error instead of returning null.
            List<PrideProxiSpectrum> spectra =
                JsonConvert.DeserializeObject<List<PrideProxiSpectrum>>(content, JsonSettings) ?? new List<PrideProxiSpectrum>();
            if (spectra.Count == 0 || spectra[0] == null)
                throw new HttpRequestException($"PRIDE PROXI returned no spectrum for USI '{usi}'.");

            return spectra[0];
        }

        /// <summary>
        /// Fetches a spectrum by its USI and returns it as an mzLib <see cref="MzSpectrum"/> — the simple case,
        /// for callers that only want the peaks. The returned object is the <see cref="PrideProxiSpectrum"/>
        /// itself (which derives from <see cref="MzSpectrum"/>), narrowed to the base type; call
        /// <see cref="GetProxiSpectrumAsync"/> instead to keep its PROXI metadata attributes in view, or
        /// <see cref="PrideArchiveExtensions.ToMsDataScan"/> to read those attributes into a full
        /// <see cref="MsDataScan"/>.
        /// </summary>
        /// <param name="usi">The Universal Spectrum Identifier.</param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns>The spectrum's peaks as an <see cref="MzSpectrum"/> (m/z ascending). Never null.</returns>
        /// <exception cref="ArgumentException">The USI is null, empty, or whitespace.</exception>
        /// <exception cref="HttpRequestException">The API returned a non-success status or an empty result for the USI.</exception>
        /// <exception cref="MzLibUtil.MzLibException">The returned spectrum's peak arrays are not parallel.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<MzSpectrum> GetSpectrumAsync(string usi, CancellationToken cancellationToken = default) =>
            await GetProxiSpectrumAsync(usi, cancellationToken).ConfigureAwait(false);

        /// <summary>
        /// Drops null elements from every collection on a deserialized <see cref="PrideProject"/>, so the
        /// type's documented "never null" guarantee covers the elements and not merely the collections.
        /// </summary>
        private static void RemoveNullElements(PrideProject project)
        {
            project.ProjectTags.RemoveAll(x => x == null);
            project.Keywords.RemoveAll(x => x == null);
            project.Countries.RemoveAll(x => x == null);
            project.Submitters.RemoveAll(x => x == null);
            project.LabPIs.RemoveAll(x => x == null);
            project.References.RemoveAll(x => x == null);
            project.Instruments.RemoveAll(x => x == null);
            project.Softwares.RemoveAll(x => x == null);
            project.ExperimentTypes.RemoveAll(x => x == null);
            project.QuantificationMethods.RemoveAll(x => x == null);
            project.Organisms.RemoveAll(x => x == null);
            project.OrganismParts.RemoveAll(x => x == null);
            project.Diseases.RemoveAll(x => x == null);
            project.IdentifiedPTMStrings.RemoveAll(x => x == null);
            project.AdditionalAttributes.RemoveAll(x => x == null);
            project.SampleAttributes.RemoveAll(x => x == null);

            // A surviving sample attribute can still hold nulls in its own value list.
            // Key and Value themselves need no guard here: an explicit "key": null or "value": null is
            // dropped by NullValueHandling.Ignore (JsonSettings), leaving the DTO's `= new()` defaults
            // standing, so neither can be null by the time this runs. That is load-bearing rather than
            // incidental -- callers dereference attribute.Key.Accession directly -- so it is pinned by
            // TryGetProjectAsync_ExplicitJsonNullKeyOrValue_DoNotClobberSampleAttributeDefaults rather
            // than re-checked here, where the null branch would be unreachable and untestable.
            foreach (PrideSampleAttribute attribute in project.SampleAttributes)
                attribute.Value.RemoveAll(x => x == null);
        }

        /// <summary>
        /// Fetches every page of a PRIDE endpoint that answers with a bare JSON array plus a
        /// <c>total_records</c> response header, and concatenates them in server order.
        /// </summary>
        /// <remarks>
        /// More than one PRIDE endpoint uses this envelope, and the termination rules below are subtle
        /// enough — and have been got wrong often enough — that a second copy of them would only be a
        /// second place for the same truncation bug to live. The rules are documented at each step
        /// rather than here, because each one exists to defend against a specific observed behavior.
        /// </remarks>
        /// <typeparam name="T">The element type of a single page.</typeparam>
        /// <param name="requestUriForPage">Builds the request URI for a zero-based page index.</param>
        /// <param name="pageSize">
        /// The page size that was requested. Used only by the fallback that runs when the response
        /// carries no usable <c>total_records</c> header.
        /// </param>
        /// <param name="subject">
        /// Names what is being paged, for error messages — e.g. <c>accession 'PXD012345'</c>. It is
        /// interpolated as written, so each endpoint's failures stay as specific as they were when this
        /// loop lived inside that endpoint's own method.
        /// </param>
        /// <param name="cancellationToken">Cancels the (possibly multi-page) fetch.</param>
        /// <returns>Every element the endpoint served, across all pages. Never null.</returns>
        private async Task<List<T>> GetAllPagesAsync<T>(Func<int, string> requestUriForPage, int pageSize,
            string subject, CancellationToken cancellationToken) where T : class
        {
            var items = new List<T>();
            string previousPageIdentity = null;
            int page = 0;
            int largestPageSeen = 0;
            bool verifyingTail = false;

            while (true)
            {
                cancellationToken.ThrowIfCancellationRequested();
                string requestUri = requestUriForPage(page);

                // Only this page's request is retried, never the fetch: restarting would re-read pages
                // already held and widen the live-index drift window SearchProjectsAsync's dedup absorbs.
                // About 6% of search requests stall for ~93 s, so without this one slow page lost the lot.
                using HttpResponseMessage response = await GetSuccessAsync(requestUri,
                    r => $"PRIDE Archive request failed with status {(int)r.StatusCode} {r.ReasonPhrase} for '{requestUri}'.",
                    cancellationToken).ConfigureAwait(false);

                string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);
                List<T> pageItems =
                    JsonConvert.DeserializeObject<List<T>>(content, JsonSettings) ?? new List<T>();

                // A null ELEMENT ("[null]") deserializes to a null entry that carries no record at all.
                // It is dropped below, mirroring what RemoveNullElements does for PrideProject, so
                // neither this loop nor a caller dereferences it. Nothing is lost: a null is not a
                // record. The count of what the SERVER actually served is taken first, because that --
                // not what survives the tail-probe dedup further down -- is what says how full a page is.
                int servedCount = pageItems.Count(x => x is not null);

                // An empty page ends the fetch. This is also the backstop for a total_records that
                // overstates what the server will actually serve. Paging past the end returns an EMPTY
                // JSON ARRAY -- "[]" from the file manifest, "[ ]" from search -- with the total_records
                // header still present (re-verified live 2026-08-21; an earlier note here claimed a
                // zero-byte body with no header, which was wrong for both endpoints. A zero-byte form
                // does exist, but only far past the end, beyond roughly offset 40 000, which no current
                // result set is large enough to reach). Every form deserializes to an empty list, the
                // zero-byte one via the null-coalesce above. An
                // overstatement is therefore accepted as the server's own correction rather than
                // reported as a shortfall -- there is no way to distinguish it from a result set whose
                // size changed mid-fetch, and throwing would fail a caller who asked for
                // nothing unreasonable.
                if (servedCount == 0)
                    break; // no (more) records

                // Every entry the server returns is kept. A page may legitimately repeat a value that
                // looks like a key (a manifest may list the same leaf name under different
                // publicFileLocations, and an entry with no fileName at all is a case DownloadFileAsync
                // already expects and guards). So paging progress is a property of the PAGE, not of the
                // individual records: deciding membership by the uniqueness of any one field would
                // silently drop real records, which is precisely the failure this loop exists to prevent.
                //
                // The page's identity is therefore its RAW RESPONSE BODY, not a projection of it.
                // A record field is the one thing this loop has already established it cannot treat as
                // an identity: a value may legitimately repeat, or be absent entirely, so two
                // consecutive pages of genuinely DIFFERENT records can share a field sequence and be
                // misread as a re-served page -- failing a correct server on exactly the result sets the
                // comment above promises to support. A separator alone does not close that gap; the
                // fields simply are not the record. The body distinguishes any two pages the server
                // actually paged, needs no per-record key the DTO does not carry, and degrades safely:
                // a server that varies its body (an embedded timestamp) merely stops this guard firing
                // and falls through to the MaxPages backstop, as it did before the guard existed.
                string previousPageBody = previousPageIdentity;
                bool pageAdvanced = content != previousPageBody;
                previousPageIdentity = content;

                // A tail probe (see the total_records branch below) is speculative: total_records has
                // already said the fetch is done, and this page is only being read in case it lied
                // downward. So it must never be able to make things worse than trusting the header
                // would have been, and there are two ways it could.
                //
                // (1) It could bring nothing new, in which case the header was right after all and the
                // fetch ends -- deliberately NOT the identical-page throw, which is reserved for a
                // server contradicting a total that says more is still to come. A server that ignores
                // `page` therefore still returns its records rather than throwing, as it did before this
                // guard existed. A byte-identical body is that case outright, and is taken here without
                // parsing anything; a body that repeats the same records in different bytes is caught by
                // the record comparison below, which sees past whitespace and property order.
                //
                // One case is NOT recovered, and the claim is limited to match: a server that both
                // ignores `page` AND varies the CONTENT of the records it re-serves (a per-record
                // timestamp, say) produces pages that are new by every measure available here, so the
                // fetch runs to the MaxPages backstop and throws where the total check used to stop it.
                // That server already defeated the identical-page guard above before this probe existed;
                // what is new is that the total check no longer rescues it. Nothing short of a record
                // key -- which this loop has established the DTO does not have -- would tell those pages
                // apart, and MaxPages remains settable by a caller who meets such a server.
                //
                // Correction, from live evidence on 2026-08-21. A body-varying server is NOT
                // hypothetical, and the paragraph above understated what it costs. On search/projects
                // PRIDE serialises the dynamic `highlights` map from an unordered hash map, so two
                // identical requests routinely differ in bytes while carrying the very same records --
                // which silently disables the identical-page guard below for that endpoint. And the
                // outcome it degrades to is NOT the MaxPages throw: a page-ignoring server whose bytes
                // vary gets its records appended a second time, items.Count reaches total, and the
                // fetch ENDS EARLY holding duplicates that look like a complete answer. SearchProjectsAsync
                // therefore deduplicates its own results on accession, which it can do because a search
                // hit carries the record key this loop does not have. The file manifest has no such key
                // and stays exposed to this -- the honest state of it, rather than a claim otherwise.
                //
                // (2) It could bring records this fetch ALREADY HOLDS. A result set that grows between
                // two requests shifts its own paging: page 0 of a 2-record set is [a, b], and page 1 of
                // the 3-record set it became is [b, c]. That is a correct server, and the comments above
                // and below both promise to tolerate it -- but appending the probe's answer wholesale
                // would hand the caller a duplicated record, which is a worse failure than the truncation
                // this probe exists to prevent, and one the caller cannot see. So a probe's records are
                // matched against the page before them and the repeats are dropped.
                //
                // Matching is on the WHOLE JSON record, not on any field of it. The distinction matters:
                // the comment above refuses to treat a field as a record's identity because a value may
                // repeat or be absent, so two different records can share one. A complete record cannot
                // be confused with a different record that way. Two byte-equal records in one manifest
                // are indistinguishable from each other anyway, so dropping one of a straddling pair
                // costs nothing a caller could act on.
                //
                // Stated plainly, because it is the one asymmetry this dedup introduces: a genuinely
                // REPEATED record survives or not depending on where the page boundary happens to fall.
                // Two byte-equal records inside one page are both kept; the same pair split across the
                // probe boundary loses one. Nothing here can tell that pair from a re-served record,
                // and the boundary is the server's to choose, so the count of an exactly-duplicated
                // record is not something a caller should read anything into.
                if (verifyingTail)
                {
                    if (!pageAdvanced)
                        break;

                    NullOutRecordsRepeatedFromPreviousPage(pageItems, content, previousPageBody);
                    if (pageItems.All(x => x is null))
                        break; // the probe held nothing this fetch did not already have
                }

                pageItems.RemoveAll(x => x is null);
                items.AddRange(pageItems);
                if (servedCount > largestPageSeen)
                    largestPageSeen = servedCount;

                if (TryGetTotalRecords(response, out long total))
                {
                    // Checked BEFORE the total comparison: a server re-serving the same page forever
                    // would otherwise push items.Count past total with duplicates and break out as if
                    // it had succeeded, hiding the fault behind a plausible-looking result.
                    //
                    // This is deliberately the strict signal -- a page byte-identical to its
                    // predecessor -- rather than "this page added nothing new". A result set that
                    // changes mid-fetch can legitimately return a page that merely overlaps the previous
                    // one, and the empty-page comment above declines to throw on exactly that
                    // ambiguity; failing only on an exactly-repeated page keeps the two consistent.
                    if (!pageAdvanced)
                        throw new MzLibException(
                            $"PRIDE Archive re-served an identical page {page} for {subject} " +
                            $"while reporting {total} total records. The server may be ignoring the page " +
                            $"parameter, or the result set may have changed mid-fetch.");

                    // The server's own record count is authoritative for stopping -- but only upward.
                    // An OVERSTATED total is already handled (the empty page above accepts it as the
                    // server's own correction). An UNDERSTATED one is the mirror failure, and it
                    // truncates the tail exactly as the capped-pageSize bug did: stop at Count >= total
                    // and every record past the reported total is dropped with no error.
                    //
                    // It is only distinguishable by asking for one more page, and only worth asking
                    // when the answer is in doubt. Truncation needs the last page to have been FULL:
                    // a server with more to give fills the page it is giving. So a page that is not
                    // full is the end of the data, the header agrees, and that stop is trusted outright.
                    // A full one is the ambiguous case -- it looks identical whether the total is honest
                    // or one page short -- and there one speculative page is fetched before believing
                    // the header.
                    //
                    // "Full" is measured against different evidence on the first page than on later
                    // ones, because the two pages carry different information.
                    //
                    // On a LATER page, the size of the largest page already served is what a full page
                    // looks like from here, and it is the only sound measure: PRIDE caps pageSize
                    // server-side and pages by the capped size, so a page of 100 against a requested 500
                    // is full, and #1102 established that the requested size cannot tell the two apart.
                    //
                    // On page 0 there is no such evidence -- the first page is trivially the largest
                    // seen, so that test is vacuously true and would probe every single-page fetch,
                    // doubling the request count of the common case. The requested pageSize is the only
                    // evidence there is, so it is used: a first page shorter than what was asked for is
                    // the server running out. That is exact whenever pageSize is at or below the server
                    // cap, which is every default call. It leaves one gap, and only one: a single-page
                    // fetch that asked for MORE than the cap gets back a capped -- and therefore
                    // possibly full -- page that looks short, so an understated total would still
                    // truncate it. That is the same above-the-cap caveat the public docstring already
                    // carries, and it buys back a request on every ordinary call.
                    //
                    // Getting this judgement wrong costs at most one request in one direction and, in
                    // the gap above, the records the header disclaimed -- never a record the header
                    // acknowledged. The probe itself cannot cost correctness at all: the guard above
                    // drops anything it repeats, so it can only ever add.
                    if (items.Count >= total)
                    {
                        bool pageCouldBeFull = page == 0
                            ? servedCount >= pageSize
                            : servedCount >= largestPageSeen;
                        if (!pageCouldBeFull)
                            break;
                        verifyingTail = true;
                    }
                    else
                    {
                        // More remain, so keep paging even though the page may look short. PRIDE caps
                        // pageSize server-side (100 as of 2026-07-23) and then pages by the capped size,
                        // so requesting 500 yields a 100-record "short" page that still has successors.
                        // Treating that as the last page silently truncated the result.
                        verifyingTail = false;
                    }
                }
                else if (servedCount < pageSize)
                {
                    // No total to trust: a short page is the last page.
                    break;
                }

                page++;
                if (page >= MaxPages)
                    throw new HttpRequestException(
                        $"PRIDE Archive paging exceeded {MaxPages} pages for {subject}; the server may be ignoring paging parameters.");
            }

            return items;
        }

        /// <summary>
        /// Nulls out the entries of a speculative tail page that the page before it already delivered,
        /// so the null sweep in <see cref="GetAllPagesAsync{T}"/> removes them before they are appended.
        /// Records are matched whole, as parsed JSON, which ignores whitespace and property order but
        /// nothing that carries meaning.
        /// </summary>
        /// <remarks>
        /// Nulling in place rather than filtering keeps this index-aligned with the raw array — the
        /// caller has not yet swept its own nulls, so element <c>i</c> here is element <c>i</c> there —
        /// and reuses the null removal that already exists instead of adding a second one.
        /// <para>
        /// Both bodies are known to be JSON arrays of the same length as their deserialized lists: the
        /// caller reached this line only by deserializing <paramref name="pageBody"/> into
        /// <paramref name="pageItems"/> element-for-element, and only on a page that follows another.
        /// Re-checking either fact here would be an unreachable branch, so the invariant is stated
        /// rather than guarded.
        /// </para>
        /// </remarks>
        private static void NullOutRecordsRepeatedFromPreviousPage<T>(List<T> pageItems, string pageBody,
            string previousPageBody) where T : class
        {
            JArray currentRecords = ParseJsonArray(pageBody);
            JArray previousRecords = ParseJsonArray(previousPageBody);

            for (int i = 0; i < pageItems.Count; i++)
            {
                if (previousRecords.Any(earlier => JToken.DeepEquals(currentRecords[i], earlier)))
                    pageItems[i] = null;
            }
        }

        /// <summary>
        /// Parses a JSON array without the date recognition Newtonsoft applies by default, so a record
        /// is compared as the server wrote it rather than as a round-trip through <see cref="DateTime"/>
        /// would render it.
        /// </summary>
        private static JArray ParseJsonArray(string json)
        {
            using var reader = new JsonTextReader(new StringReader(json)) { DateParseHandling = DateParseHandling.None };
            return JArray.Load(reader);
        }

        /// <summary>Reads the PRIDE "total_records" response header, if present and numeric.</summary>
        private static bool TryGetTotalRecords(HttpResponseMessage response, out long total)
        {
            total = 0;
            if (response.Headers.TryGetValues("total_records", out IEnumerable<string> values))
            {
                foreach (string value in values)
                {
                    if (long.TryParse(value, out total))
                        return true;
                }
            }
            return false;
        }

        /// <inheritdoc/>
        public void Dispose()
        {
            if (_disposed)
                return;
            if (_ownsHttpClient)
                _httpClient.Dispose();
            _disposed = true;
        }
    }
}
