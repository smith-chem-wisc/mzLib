using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
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
        /// How long <see cref="DownloadFileAsync"/> will wait for the NEXT bytes of a response body before
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
        /// <exception cref="HttpRequestException">The API returned a non-success status code.</exception>
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
        /// host also serves that path over HTTPS — the same fact <see cref="DownloadFileAsync"/> relies
        /// on — so the whole walk goes over this client's reused <see cref="HttpClient"/>. Subdirectories
        /// are followed; each returned file's <see cref="PrideFtpFile.RelativePath"/> is relative to the
        /// project root. Sizes are PRIDE's rounded index sizes — see <see cref="PrideFtpFile.ApproximateSizeBytes"/>.
        /// </remarks>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="cancellationToken">Cancels the (multi-request) walk.</param>
        /// <returns>Every file under the project's FTP root, subdirectories included. Never null.</returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="MzLibException">No project has that accession, or it carries no publication date to locate its FTP directory.</exception>
        /// <exception cref="HttpRequestException">The project or a directory could not be fetched (non-success status), or the directory nesting exceeded the cyclic-listing guard.</exception>
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
            await CollectFtpFilesAsync(new Uri(rootUrl), relativePrefix: "", depth: 0, files, cancellationToken)
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
        private async Task CollectFtpFilesAsync(Uri directoryUri, string relativePrefix, int depth,
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

            using HttpResponseMessage response = await _httpClient.GetAsync(directoryUri, cancellationToken).ConfigureAwait(false);
            if (!response.IsSuccessStatusCode)
                throw new HttpRequestException(
                    $"PRIDE FTP directory listing failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{directoryUri}'.");

            string html = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

            foreach ((string rawHref, string sizeText) in ParseAutoIndexRows(html))
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
                    await CollectFtpFilesAsync(childUri, relativePrefix + name, depth + 1, files, cancellationToken)
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
            using HttpResponseMessage response = await _httpClient.GetAsync(requestUri, cancellationToken).ConfigureAwait(false);

            // 404 is the ONE expected "no" and is reported as a value. It is checked before the general
            // status guard so that every other failure still throws.
            if (response.StatusCode == HttpStatusCode.NotFound)
                return (false, null);

            if (!response.IsSuccessStatusCode)
                throw new HttpRequestException(
                    $"PRIDE Archive request failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{requestUri}'.");

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
        /// <exception cref="HttpRequestException">The API returned a non-success status code.</exception>
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
        /// <exception cref="HttpRequestException">The download returned a non-success status code, or the response body delivered nothing for <see cref="BodyStallTimeout"/>.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled. A stall is NOT reported this way — see <see cref="BodyStallTimeout"/>.</exception>
        public async Task<string> DownloadFileAsync(PrideArchiveFile file, string destinationDirectory,
            bool overwrite = true, CancellationToken cancellationToken = default)
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

            string destinationPath = Path.Combine(destinationDirectory, safeFileName);

            // Cheap resume: an already-present destination is left untouched. This runs before URL
            // resolution and directory creation so skipping a downloaded file never fails on a file
            // that has no HTTPS location (e.g. Aspera-only) or does needless filesystem work.
            if (!overwrite && File.Exists(destinationPath))
                return destinationPath;

            string url = file.GetHttpsDownloadUrl(); // throws NotSupportedException if unreachable over HTTPS

            Directory.CreateDirectory(destinationDirectory);

            using HttpResponseMessage response =
                await _httpClient.GetAsync(url, HttpCompletionOption.ResponseHeadersRead, cancellationToken).ConfigureAwait(false);
            if (!response.IsSuccessStatusCode)
                throw new HttpRequestException(
                    $"PRIDE download failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{url}'.");

            string partialPath = destinationPath + ".partial";
            try
            {
                using (var fileStream = new FileStream(partialPath, FileMode.Create, FileAccess.Write, FileShare.None))
                using (Stream httpStream = await response.Content.ReadAsStreamAsync(cancellationToken).ConfigureAwait(false))
                {
                    // Not Stream.CopyToAsync: it would inherit the very absence of a read deadline that
                    // BodyStallTimeout exists to supply, which is how the body escaped every timeout here.
                    await CopyUntilStalledAsync(httpStream, fileStream, url, cancellationToken).ConfigureAwait(false);
                }
                File.Move(partialPath, destinationPath, overwrite: true);
            }
            finally
            {
                if (File.Exists(partialPath))
                    File.Delete(partialPath);
            }

            return destinationPath;
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
        private async Task CopyUntilStalledAsync(Stream source, Stream destination, string url,
            CancellationToken cancellationToken)
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
                            $"The PRIDE response body for '{url}' delivered nothing for {BodyStallTimeout}.", e);
                    }
                }

                if (read == 0)
                    return;

                await destination.WriteAsync(buffer.AsMemory(0, read), cancellationToken).ConfigureAwait(false);
            }
        }

        /// <summary>
        /// Downloads a project's files to <paramref name="destinationDirectory"/>, optionally filtered. This is
        /// the convenience over <see cref="GetProjectFilesAsync"/> + <see cref="DownloadFileAsync"/>: it fetches
        /// the manifest, applies <paramref name="filter"/>, and downloads each selected file in turn.
        /// </summary>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="destinationDirectory">The directory to write into; created if it does not exist.</param>
        /// <param name="filter">
        /// An optional predicate selecting which files to download (e.g. by category or extension — see
        /// <see cref="PrideArchiveExtensions"/>). When null, every file in the manifest is downloaded.
        /// </param>
        /// <param name="overwrite">Passed through to <see cref="DownloadFileAsync"/>; default true.</param>
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
            using HttpResponseMessage response = await _httpClient.GetAsync(requestUri, cancellationToken).ConfigureAwait(false);

            if (!response.IsSuccessStatusCode)
                throw new HttpRequestException(
                    $"PRIDE PROXI request failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{requestUri}'.");

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
                using HttpResponseMessage response = await _httpClient.GetAsync(requestUri, cancellationToken).ConfigureAwait(false);

                if (!response.IsSuccessStatusCode)
                    throw new HttpRequestException(
                        $"PRIDE Archive request failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{requestUri}'.");

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
