using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Threading;
using System.Threading.Tasks;
using Newtonsoft.Json;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A client for the PRIDE Archive REST API (https://www.ebi.ac.uk/pride/ws/archive/v3/) that makes
    /// EBI PRIDE proteomics dataset data available to mzLib. The initial capability, given a project
    /// accession (e.g. "PXD012345"), is to return the complete manifest of that project's files.
    /// </summary>
    /// <remarks>
    /// Holds one reusable <see cref="HttpClient"/> and is <see cref="IDisposable"/>. Use one instance
    /// per unit of work and dispose it. A constructor overload accepts an <see cref="HttpClient"/> to
    /// support testing and custom configuration.
    /// </remarks>
    public sealed class PrideArchiveClient : IDisposable
    {
        /// <summary>The base address of the PRIDE Archive REST API (v3).</summary>
        public const string DefaultBaseAddress = "https://www.ebi.ac.uk/pride/ws/archive/v3/";

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
            if (safeFileName != file.FileName)
                throw new ArgumentException(
                    $"The PRIDE file name '{file.FileName}' is not a bare file name; refusing to write outside the destination directory.",
                    nameof(file));

            string url = file.GetHttpsDownloadUrl(); // throws NotSupportedException if unreachable over HTTPS

            Directory.CreateDirectory(destinationDirectory);
            string destinationPath = Path.Combine(destinationDirectory, safeFileName);

            if (!overwrite && File.Exists(destinationPath))
                return destinationPath; // already downloaded; leave it untouched

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
