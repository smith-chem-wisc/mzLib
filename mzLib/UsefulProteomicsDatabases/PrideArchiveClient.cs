using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
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
        /// Returns the complete manifest of files belonging to a PRIDE Archive project. All pages are
        /// fetched and concatenated; paging is an implementation detail hidden from the caller.
        /// </summary>
        /// <param name="accession">The PRIDE project accession, e.g. "PXD012345".</param>
        /// <param name="pageSize">Files requested per page (default 100). The full manifest is returned regardless.</param>
        /// <param name="cancellationToken">Cancels the (possibly multi-page) fetch.</param>
        /// <returns>
        /// The project's files. Empty if the project has no files or the accession is unknown (PRIDE
        /// returns an empty result for an unknown accession). Never null.
        /// </returns>
        /// <exception cref="ArgumentException">The accession is null, empty, or whitespace.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The page size is not positive.</exception>
        /// <exception cref="HttpRequestException">The API returned a non-success status code.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<List<PrideArchiveFile>> GetProjectFilesAsync(string accession, int pageSize = 100,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A PRIDE project accession is required.", nameof(accession));
            if (pageSize <= 0)
                throw new ArgumentOutOfRangeException(nameof(pageSize), "Page size must be positive.");

            var files = new List<PrideArchiveFile>();
            int page = 0;

            while (true)
            {
                cancellationToken.ThrowIfCancellationRequested();
                string requestUri = $"projects/{Uri.EscapeDataString(accession)}/files?pageSize={pageSize}&page={page}";
                using HttpResponseMessage response = await _httpClient.GetAsync(requestUri, cancellationToken).ConfigureAwait(false);

                if (!response.IsSuccessStatusCode)
                    throw new HttpRequestException(
                        $"PRIDE Archive request failed with status {(int)response.StatusCode} {response.ReasonPhrase} for '{requestUri}'.");

                string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);
                List<PrideArchiveFile> pageFiles =
                    JsonConvert.DeserializeObject<List<PrideArchiveFile>>(content, JsonSettings) ?? new List<PrideArchiveFile>();

                if (pageFiles.Count == 0)
                    break; // no (more) files

                files.AddRange(pageFiles);

                // Stop when we have collected everything the server reported.
                if (TryGetTotalRecords(response, out long total) && files.Count >= total)
                    break;
                // Safety net: a short page is the last page even if the total header is missing.
                if (pageFiles.Count < pageSize)
                    break;

                page++;
                if (page >= MaxPages)
                    throw new HttpRequestException(
                        $"PRIDE Archive paging exceeded {MaxPages} pages for accession '{accession}'; the server may be ignoring paging parameters.");
            }

            return files;
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
            foreach (PrideSampleAttribute attribute in project.SampleAttributes)
                attribute.Value.RemoveAll(x => x == null);
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
