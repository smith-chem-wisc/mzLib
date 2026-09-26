using System;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A paper as Europe PMC knows it: its identifiers and whether its full text and supplements can be fetched.
    /// Strings are never null.
    /// </summary>
    public sealed class EuropePmcArticle
    {
        /// <summary>The PubMed id, or empty.</summary>
        public string PubMedId { get; init; } = string.Empty;

        /// <summary>The PubMed Central id (<c>PMC6910996</c>), or empty when the paper is not in PMC.</summary>
        public string PmcId { get; init; } = string.Empty;

        /// <summary>The DOI, or empty.</summary>
        public string Doi { get; init; } = string.Empty;

        /// <summary>The title, or empty.</summary>
        public string Title { get; init; } = string.Empty;

        /// <summary>True when Europe PMC marks the paper open access.</summary>
        public bool IsOpenAccess { get; init; }

        /// <summary>
        /// True when Europe PMC holds the paper's full text, so <see cref="EuropePmcClient.TryDownloadFullTextXmlAsync"/>
        /// and <see cref="EuropePmcClient.TryDownloadSupplementaryFilesAsync"/> may ask for it.
        /// </summary>
        public bool HasFullText { get; init; }
    }

    /// <summary>
    /// A client for the Europe PMC REST API (https://www.ebi.ac.uk/europepmc/webservices/rest/): resolves a PRIDE
    /// project's reference (PubMed id or DOI) to the paper, and downloads its full text (JATS XML) and supplementary
    /// files. Built for drafting SDRFs from publications (the sdrf project's design SAMPLE-EVIDENCE.md, E4); its
    /// caller is the drafting pipeline, outside mzLib, so it is public.
    /// </summary>
    /// <remarks>
    /// Follows <see cref="PrideArchiveClient"/>. An absence is a value: no such paper, no open full text, no
    /// supplements are <c>Found = false</c>. The service being unreachable, timing out, rate-limiting (408/429) or
    /// answering 5xx is an <see cref="HttpRequestException"/>, which <c>ExternalServiceTestHelper</c> reads as an
    /// outage; any other refusal (401/403/400) is an <see cref="MzLibException"/>, as in <see cref="ProteinDbRetriever"/>,
    /// so a changed contract fails instead of skipping. Messages name the endpoint, never the URL.
    /// <para>
    /// <b>Two Europe PMC habits it guards against.</b> A closed article's full text is answered with HTTP 500 -- an
    /// outage by every rule above -- so full text is only requested when <see cref="EuropePmcArticle.HasFullText"/>
    /// says it exists. And "no supplements" is a 200 carrying an XML error page, so a supplement download is kept
    /// only if it is a zip archive; an error body is never written to disk.
    /// </para>
    /// <para>
    /// A downloaded file already on disk is reused unless <c>overwrite</c> is true (the PRIDE client's cheap resume),
    /// which makes a folder of downloads a cache. Downloads stream to <c>.partial</c> and are moved into place.
    /// </para>
    /// </remarks>
    public sealed class EuropePmcClient : IDisposable
    {
        /// <summary>The base address of the Europe PMC REST API.</summary>
        public const string DefaultBaseAddress = "https://www.ebi.ac.uk/europepmc/webservices/rest/";

        private static readonly JsonSerializerSettings JsonSettings = new() { NullValueHandling = NullValueHandling.Ignore };

        private readonly HttpClient _httpClient;
        private readonly bool _ownsHttpClient;
        private bool _disposed;

        /// <summary>How long a download waits for the NEXT bytes before abandoning the transfer (see <see cref="PrideArchiveClient.BodyStallTimeout"/>).</summary>
        public TimeSpan BodyStallTimeout { get; init; } = TimeSpan.FromMinutes(2);

        /// <summary>Creates a client with its own <see cref="HttpClient"/>.</summary>
        public EuropePmcClient()
            : this(new HttpClient { BaseAddress = new Uri(DefaultBaseAddress), Timeout = TimeSpan.FromSeconds(100) }, ownsHttpClient: true)
        {
        }

        /// <summary>
        /// Creates a client over a supplied <see cref="HttpClient"/>; the caller keeps ownership. A client with no
        /// <see cref="HttpClient.BaseAddress"/> gets <see cref="DefaultBaseAddress"/>.
        /// </summary>
        public EuropePmcClient(HttpClient httpClient)
            : this(httpClient, ownsHttpClient: false)
        {
        }

        private EuropePmcClient(HttpClient httpClient, bool ownsHttpClient)
        {
            _httpClient = httpClient ?? throw new ArgumentNullException(nameof(httpClient));
            _ownsHttpClient = ownsHttpClient;
            if (_httpClient.BaseAddress == null)
                _httpClient.BaseAddress = new Uri(DefaultBaseAddress);
        }

        /// <summary>Finds the paper a PRIDE project cites, by PubMed id first, else DOI.</summary>
        /// <exception cref="ArgumentNullException"><paramref name="reference"/> is null.</exception>
        /// <exception cref="ArgumentException">The reference has neither a PubMed id nor a DOI.</exception>
        public Task<(bool Found, EuropePmcArticle Article)> TryFindArticleAsync(PrideReference reference, CancellationToken cancellationToken = default)
        {
            if (reference == null) throw new ArgumentNullException(nameof(reference));
            return TryFindArticleAsync(reference.PubmedId, reference.Doi, cancellationToken);
        }

        /// <summary>
        /// Finds a paper by PubMed id (a positive number) or, when there is none, by DOI. <c>Found</c> is false only
        /// when Europe PMC answered and has no such paper.
        /// </summary>
        /// <exception cref="ArgumentException">Neither a PubMed id nor a DOI was given.</exception>
        /// <exception cref="HttpRequestException">Europe PMC was unreachable, timed out, rate-limited or answered 5xx.</exception>
        /// <exception cref="MzLibException">Europe PMC refused the request otherwise, or answered with something that is not a search result.</exception>
        public async Task<(bool Found, EuropePmcArticle Article)> TryFindArticleAsync(int pubMedId, string doi, CancellationToken cancellationToken = default)
        {
            string query = pubMedId > 0 ? $"EXT_ID:{pubMedId} AND SRC:MED"
                : !string.IsNullOrWhiteSpace(doi) ? $"DOI:\"{doi.Trim()}\""
                : throw new ArgumentException("A PubMed id or a DOI is required.", nameof(doi));
            cancellationToken.ThrowIfCancellationRequested();

            string requestUri = "search?format=json&resultType=lite&query=" + Uri.EscapeDataString(query);
            using HttpResponseMessage response = await _httpClient.GetAsync(requestUri, cancellationToken).ConfigureAwait(false);
            ThrowIfFailed(response, "search");
            string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

            JObject json;
            try { json = JObject.Parse(content); }
            catch (JsonReaderException e) { throw new MzLibException("Europe PMC search answered with something that is not JSON.", e); }
            if (json["resultList"]?["result"] is not JArray results)
                throw new MzLibException("Europe PMC search answered without a result list.");
            var hit = results.OfType<JObject>().FirstOrDefault();
            if (hit == null) return (false, null);

            string pmcid = (string)hit["pmcid"] ?? string.Empty;
            bool inEpmc = string.Equals((string)hit["inEPMC"], "Y", StringComparison.OrdinalIgnoreCase);
            bool listed = hit["fullTextIdList"]?["fullTextId"] is JArray ids && ids.Any(i => string.Equals((string)i, pmcid, StringComparison.OrdinalIgnoreCase));
            return (true, new EuropePmcArticle
            {
                PubMedId = (string)hit["pmid"] ?? string.Empty,
                PmcId = pmcid,
                Doi = (string)hit["doi"] ?? string.Empty,
                Title = (string)hit["title"] ?? string.Empty,
                IsOpenAccess = string.Equals((string)hit["isOpenAccess"], "Y", StringComparison.OrdinalIgnoreCase),
                HasFullText = pmcid.Length > 0 && (inEpmc || listed)
            });
        }

        /// <summary>
        /// Saves the paper's full text as <c>&lt;PMCID&gt;.xml</c> (JATS) in <paramref name="directory"/>. <c>Found</c> is
        /// false, and nothing is requested, when the article has no full text in Europe PMC.
        /// </summary>
        /// <exception cref="MzLibException">Europe PMC answered 200 with something that is not XML.</exception>
        public async Task<(bool Found, string Path)> TryDownloadFullTextXmlAsync(EuropePmcArticle article, string directory,
            bool overwrite = false, CancellationToken cancellationToken = default)
        {
            if (!Usable(article, directory)) return (false, null);
            string path = System.IO.Path.Combine(directory, article.PmcId + ".xml");
            if (!overwrite && File.Exists(path)) return (true, path);

            byte[] body = await DownloadAsync($"{Uri.EscapeDataString(article.PmcId)}/fullTextXML", "full text", cancellationToken).ConfigureAwait(false);
            if (body == null) return (false, null);
            if (!LooksLikeXml(body))
                throw new MzLibException($"Europe PMC full text for {article.PmcId} is not XML.");
            Save(path, body);
            return (true, path);
        }

        /// <summary>
        /// Saves the paper's supplementary files as one zip, <c>&lt;PMCID&gt;_supplementary.zip</c>, in
        /// <paramref name="directory"/>. <c>Found</c> is false when the article has none in Europe PMC -- which Europe
        /// PMC reports as a 200 carrying an error page, and which is never written to disk.
        /// </summary>
        public async Task<(bool Found, string Path)> TryDownloadSupplementaryFilesAsync(EuropePmcArticle article, string directory,
            bool overwrite = false, CancellationToken cancellationToken = default)
        {
            if (!Usable(article, directory)) return (false, null);
            string path = System.IO.Path.Combine(directory, article.PmcId + "_supplementary.zip");
            if (!overwrite && File.Exists(path)) return (true, path);

            byte[] body = await DownloadAsync($"{Uri.EscapeDataString(article.PmcId)}/supplementaryFiles?includeInlineImage=false",
                "supplementary files", cancellationToken).ConfigureAwait(false);
            // A zip starts "PK"; anything else (the XML errorBean, an HTML page) means there are none.
            if (body == null || body.Length < 4 || body[0] != (byte)'P' || body[1] != (byte)'K') return (false, null);
            Save(path, body);
            return (true, path);
        }

        private static bool Usable(EuropePmcArticle article, string directory)
        {
            if (article == null) throw new ArgumentNullException(nameof(article));
            if (string.IsNullOrWhiteSpace(directory)) throw new ArgumentException("A destination directory is required.", nameof(directory));
            // A closed article's full text is answered with HTTP 500: never ask for it (see the remarks).
            if (!article.HasFullText || string.IsNullOrWhiteSpace(article.PmcId)) return false;
            if (article.PmcId.IndexOfAny(System.IO.Path.GetInvalidFileNameChars()) >= 0 || article.PmcId.Contains(".."))
                throw new ArgumentException($"'{article.PmcId}' is not a PMC id.", nameof(article));
            return true;
        }

        /// <summary>The body, or null on 404. Throws on every other failure (see the remarks).</summary>
        private async Task<byte[]> DownloadAsync(string requestUri, string what, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            using HttpResponseMessage response =
                await _httpClient.GetAsync(requestUri, HttpCompletionOption.ResponseHeadersRead, cancellationToken).ConfigureAwait(false);
            if (response.StatusCode == HttpStatusCode.NotFound) return null;
            ThrowIfFailed(response, what);
            using var buffer = new MemoryStream();
            using (Stream body = await response.Content.ReadAsStreamAsync(cancellationToken).ConfigureAwait(false))
                await CopyUntilStalledAsync(body, buffer, what, cancellationToken).ConfigureAwait(false);
            return buffer.ToArray();
        }

        private static void ThrowIfFailed(HttpResponseMessage response, string what)
        {
            if (response.IsSuccessStatusCode) return;
            int status = (int)response.StatusCode;
            string message = $"Europe PMC {what} request failed with status {status} {response.ReasonPhrase}.";
            if (status == 408 || status == 429 || status >= 500) throw new HttpRequestException(message, null, response.StatusCode);
            throw new MzLibException(message);
        }

        /// <summary>Writes to <c>.partial</c> and moves into place, so a half-written file is never taken for a download.</summary>
        private static void Save(string path, byte[] body)
        {
            Directory.CreateDirectory(System.IO.Path.GetDirectoryName(path));
            string partial = path + ".partial";
            try
            {
                File.WriteAllBytes(partial, body);
                File.Move(partial, path, overwrite: true);
            }
            finally
            {
                if (File.Exists(partial)) File.Delete(partial);
            }
        }

        private static bool LooksLikeXml(byte[] body)
        {
            int i = body.Length >= 3 && body[0] == 0xEF && body[1] == 0xBB && body[2] == 0xBF ? 3 : 0;
            while (i < body.Length && char.IsWhiteSpace((char)body[i])) i++;
            return i < body.Length && body[i] == (byte)'<';
        }

        private async Task CopyUntilStalledAsync(Stream source, Stream destination, string what, CancellationToken cancellationToken)
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
                    catch (OperationCanceledException e) when (stallWindow.IsCancellationRequested && !cancellationToken.IsCancellationRequested)
                    {
                        throw new HttpRequestException($"The Europe PMC {what} response delivered nothing for {BodyStallTimeout}.", e);
                    }
                }
                if (read == 0) return;
                await destination.WriteAsync(buffer.AsMemory(0, read), cancellationToken).ConfigureAwait(false);
            }
        }

        /// <inheritdoc/>
        public void Dispose()
        {
            if (_disposed) return;
            if (_ownsHttpClient) _httpClient.Dispose();
            _disposed = true;
        }
    }
}
