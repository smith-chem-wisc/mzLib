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
        /// <exception cref="HttpRequestException">The download returned a non-success status code.</exception>
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
                    await httpStream.CopyToAsync(fileStream, cancellationToken).ConfigureAwait(false);
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
                // Drop it here, mirroring what RemoveNullElements does for PrideProject, so neither this
                // loop nor a caller dereferences it. Nothing is lost: a null is not a record.
                pageItems.RemoveAll(x => x is null);

                // An empty page ends the fetch. This is also the backstop for a total_records that
                // overstates what the server will actually serve: paging past the end returns a
                // zero-byte body (verified live 2026-07-23), which deserializes to an empty list. An
                // overstatement is therefore accepted as the server's own correction rather than
                // reported as a shortfall -- there is no way to distinguish it from a result set whose
                // size changed mid-fetch, and throwing would fail a caller who asked for
                // nothing unreasonable.
                if (pageItems.Count == 0)
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
                bool pageAdvanced = content != previousPageIdentity;
                previousPageIdentity = content;

                items.AddRange(pageItems);

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

                    // The server's own record count is authoritative: stop once we have all of it.
                    if (items.Count >= total)
                        break;
                    // More remain, so keep paging even though the page may look short. PRIDE caps
                    // pageSize server-side (100 as of 2026-07-23) and then pages by the capped size,
                    // so requesting 500 yields a 100-record "short" page that still has successors.
                    // Treating that as the last page silently truncated the result.
                }
                else if (pageItems.Count < pageSize)
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
