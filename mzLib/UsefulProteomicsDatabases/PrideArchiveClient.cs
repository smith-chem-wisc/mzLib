using System;
using System.Collections.Generic;
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
