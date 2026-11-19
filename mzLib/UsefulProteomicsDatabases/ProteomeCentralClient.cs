using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A client for ProteomeCentral, the ProteomeXchange consortium's catalogue
    /// (https://proteomecentral.proteomexchange.org). One search reaches every ProteomeXchange repository at once:
    /// PRIDE, MassIVE, jPOST, iProX and PanoramaPublic.
    /// </summary>
    /// <remarks>
    /// <b>It knows only datasets with a PXD.</b> A MassIVE dataset that was never announced to ProteomeXchange is
    /// invisible here (48 of MassIVE's 67 "zebrafish" datasets, 2026-11-10); <see cref="MassiveClient.SearchDatasetsAsync"/>
    /// finds those.
    /// <para>
    /// The failure contract is the other repository clients': an outage is an <see cref="HttpRequestException"/>,
    /// retried first (3 times, after 5, 20 and 60 s) and then marked for <see cref="ServiceFailure.IsTransient"/>; an
    /// answer that breaks the contract is an <see cref="MzLibException"/>. Holds one <see cref="HttpClient"/>; dispose it.
    /// </para>
    /// </remarks>
    public sealed class ProteomeCentralClient : IDisposable
    {
        /// <summary>The base address of ProteomeCentral's PROXI service.</summary>
        public const string DefaultBaseAddress = "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/";

        /// <summary>
        /// The longest keyword <see cref="SearchDatasetsAsync"/> will send: PRIDE's measured limit. ProteomeCentral
        /// answered a 1,000-character search with 200 (2026-11-10); nothing longer was tried.
        /// </summary>
        public const int MaxKeywordLength = 1000;

        // The columns of a compact search answer, by the names in its result_set.datasets_title_list. Read by name, not
        // by position, so a column that moves is still read right and one that vanishes is null rather than another's value.
        private const string IdentifierColumn = "dataset identifier";
        private const string TitleColumn = "title";
        private const string RepositoryColumn = "repository";
        private const string SpeciesColumn = "species";
        private const string SdrfColumn = "SDRF";
        private const string FilesColumn = "files (raw/total)";
        private const string InstrumentColumn = "instrument";
        private const string PublicationsColumn = "publications";
        private const string AnnounceDateColumn = "announce date";
        private const string KeywordsColumn = "keywords";

        private readonly HttpClient _httpClient;
        private readonly bool _ownsHttpClient;
        private readonly DownloadEngine _engine;
        private bool _disposed;

        /// <summary>
        /// A hard upper bound on the pages <see cref="SearchDatasetsAsync"/> fetches, against a server that ignores the
        /// paging parameters; exceeding it throws.
        /// </summary>
        public int MaxPages { get; init; } = 10000;

        /// <summary>How many times a transient failure is retried after the first attempt (3). Tests set 0 to pin a single attempt.</summary>
        internal int MaxRetries { get => _engine.MaxRetries; init => _engine.MaxRetries = value; }

        /// <summary>Waits out a backoff. Tests replace it to record the schedule instead of sleeping through it.</summary>
        internal Func<TimeSpan, CancellationToken, Task> RetryDelay { get => _engine.RetryDelay; init => _engine.RetryDelay = value; }

        /// <summary>Creates a client with its own <see cref="HttpClient"/>.</summary>
        public ProteomeCentralClient()
            : this(new HttpClient { Timeout = TimeSpan.FromSeconds(100) }, ownsHttpClient: true)
        {
        }

        /// <summary>
        /// Creates a client over a supplied <see cref="HttpClient"/>. The caller retains ownership: it is NOT disposed by
        /// <see cref="Dispose"/>. Requests are absolute, so its base address is unused.
        /// </summary>
        /// <param name="httpClient">The HTTP client to use. Must not be null.</param>
        public ProteomeCentralClient(HttpClient httpClient)
            : this(httpClient, ownsHttpClient: false)
        {
        }

        private ProteomeCentralClient(HttpClient httpClient, bool ownsHttpClient)
        {
            _httpClient = httpClient ?? throw new ArgumentNullException(nameof(httpClient));
            _ownsHttpClient = ownsHttpClient;
            _engine = new DownloadEngine(_httpClient, "ProteomeCentral");
        }

        /// <summary>
        /// Searches every ProteomeXchange repository's datasets, by keyword, by filter, or both.
        /// </summary>
        /// <remarks>
        /// Each filter in <paramref name="filter"/> was measured to filter (2026-11-10): asking for a facet's top value
        /// returned exactly that facet's count. Only those parameters are ever sent, because ProteomeCentral ignores a
        /// parameter it does not know and answers with every dataset. A VALUE it does not know returns nothing, so a
        /// misspelt value gives an empty list rather than everything. Values are ProteomeCentral's display names
        /// ("Homo sapiens", "Q Exactive", "MassIVE", "curated valid SDRF").
        /// <para>
        /// Every page is fetched; <paramref name="pageSize"/> only sets how many datasets each request asks for. The
        /// search ends at a page shorter than asked for (an empty one included), so a total ProteomeCentral understates
        /// cannot cut it short. A dataset served twice as the catalogue moves is kept once, by its PXD.
        /// </para>
        /// <para>
        /// The keyword matches titles and submitter keywords, not descriptions (aging's survey, 2026-11-10). Its
        /// <c>resultType=full</c> returns the same rows, so a description or modifications need the dataset's own record:
        /// <see cref="MassiveClient.GetDatasetAsync"/> or <see cref="PrideArchiveClient.GetProjectAsync"/>. To combine
        /// this search with MassIVE's, see <see cref="RepositorySearch.Merge"/>.
        /// </para>
        /// </remarks>
        /// <param name="keyword">The text to search for; may be null or blank when <paramref name="filter"/> sets a value.</param>
        /// <param name="filter">Optional filters; every value set must match for a dataset to be returned.</param>
        /// <param name="pageSize">How many datasets each request asks for. Positive.</param>
        /// <param name="cancellationToken">Cancels the (multi-request) search.</param>
        /// <returns>Every dataset found, in ProteomeCentral's order, each once. Empty when nothing matches. Never null.</returns>
        /// <exception cref="ArgumentException">Neither a keyword nor a filter value is given; the keyword is longer than <see cref="MaxKeywordLength"/>; or a filter value is blank. No request is made.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="pageSize"/> is not positive.</exception>
        /// <exception cref="MzLibException">An answer is not JSON, has no column list or no "dataset identifier" column, has a row with no identifier, or served a full page holding only datasets already served (the paging is being ignored).</exception>
        /// <exception cref="HttpRequestException">ProteomeCentral was unreachable or answered a non-success status, after the retries for a transient one; or the search ran past <see cref="MaxPages"/>.</exception>
        /// <exception cref="OperationCanceledException">The search was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<IReadOnlyList<ProteomeCentralDataset>> SearchDatasetsAsync(string keyword,
            ProteomeCentralSearchFilter filter = null, int pageSize = 100, CancellationToken cancellationToken = default)
        {
            bool hasKeyword = !string.IsNullOrWhiteSpace(keyword);
            IReadOnlyList<(string Name, string Value)> filters = filter?.Parameters() ?? Array.Empty<(string, string)>();
            if (!hasKeyword && filters.Count == 0)
                throw new ArgumentException("A search keyword, a filter value, or both are required.", nameof(keyword));
            if (hasKeyword && keyword.Length > MaxKeywordLength)
                throw new ArgumentException(
                    $"The search keyword is {keyword.Length} characters; at most {MaxKeywordLength} are sent.", nameof(keyword));
            if (pageSize <= 0)
                throw new ArgumentOutOfRangeException(nameof(pageSize), "Page size must be positive.");
            ObjectDisposedException.ThrowIf(_disposed, this);

            var query = new StringBuilder();
            if (hasKeyword)
                query.Append("&search=").Append(Uri.EscapeDataString(keyword));
            foreach ((string name, string value) in filters)
                query.Append('&').Append(name).Append('=').Append(Uri.EscapeDataString(value));
            string subject = hasKeyword ? $"keyword '{keyword}'" : "the filter";

            var results = new List<ProteomeCentralDataset>();
            var seen = new HashSet<string>(StringComparer.Ordinal);

            for (int page = 1; ; page++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (page > MaxPages)
                    throw new HttpRequestException(
                        $"ProteomeCentral search paging exceeded {MaxPages} pages for {subject}; the server may be ignoring the paging parameters.");

                string requestUri = $"{DefaultBaseAddress}datasets?resultType=compact&pageSize={pageSize}&pageNumber={page}{query}";
                string described = $"page {page} of ProteomeCentral's search for {subject}";
                using HttpResponseMessage response = await _engine.GetSuccessAsync(requestUri,
                    r => $"ProteomeCentral search failed with status {(int)r.StatusCode} {r.ReasonPhrase} for {described}.",
                    cancellationToken).ConfigureAwait(false);
                string json = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

                List<ProteomeCentralDataset> rows = ReadPage(json, described);
                int added = 0;
                foreach (ProteomeCentralDataset row in rows)
                {
                    if (seen.Add(row.Accession))
                    {
                        results.Add(row);
                        added++;
                    }
                }

                if (rows.Count < pageSize)
                    break;
                if (page > 1 && added == 0)
                    throw new MzLibException(
                        $"ProteomeCentral served {described} with only datasets it had already served; it is ignoring the paging parameters.");
            }

            return results;
        }

        /// <summary>The datasets of one compact search answer, read by column name.</summary>
        internal static List<ProteomeCentralDataset> ReadPage(string json, string described)
        {
            JObject answer;
            try
            {
                using var reader = new JsonTextReader(new System.IO.StringReader(json)) { DateParseHandling = DateParseHandling.None };
                answer = JToken.ReadFrom(reader) as JObject;
            }
            catch (JsonReaderException e)
            {
                throw new MzLibException($"ProteomeCentral's answer to {described} is not valid JSON.", e);
            }

            if (answer?["result_set"]?["datasets_title_list"] is not JArray titles)
                throw new MzLibException($"ProteomeCentral's answer to {described} has no column list (result_set.datasets_title_list).");

            var column = new Dictionary<string, int>(StringComparer.Ordinal);
            for (int i = 0; i < titles.Count; i++)
                if (titles[i].Type == JTokenType.String)
                    column.TryAdd((string)titles[i], i);
            if (!column.ContainsKey(IdentifierColumn))
                throw new MzLibException(
                    $"ProteomeCentral's answer to {described} has no \"{IdentifierColumn}\" column, so its rows cannot be named.");

            var rows = new List<ProteomeCentralDataset>();
            if (answer["datasets"] is not JArray datasets)
                return rows; // a page past the last carries no rows

            foreach (JToken token in datasets)
            {
                if (token is not JArray row)
                    throw new MzLibException($"A row of ProteomeCentral's answer to {described} is not a list of values.");

                string Cell(string name) =>
                    column.TryGetValue(name, out int i) && i < row.Count && row[i].Type is not JTokenType.Null and not JTokenType.Undefined
                        ? WebUtility.HtmlDecode(row[i].Type == JTokenType.String ? (string)row[i] : row[i].ToString(Formatting.None))
                        : null;

                string accession = Cell(IdentifierColumn);
                if (string.IsNullOrWhiteSpace(accession))
                    throw new MzLibException($"A row of ProteomeCentral's answer to {described} has no dataset identifier.");

                rows.Add(new ProteomeCentralDataset
                {
                    Accession = accession,
                    Title = Cell(TitleColumn),
                    Repository = Cell(RepositoryColumn),
                    Species = Cell(SpeciesColumn),
                    Sdrf = Cell(SdrfColumn),
                    Files = Cell(FilesColumn),
                    Instrument = Cell(InstrumentColumn),
                    PublicationText = PublicationTextOf(Cell(PublicationsColumn)),
                    PubMedIds = LinkedIds(Cell(PublicationsColumn), PubMedLink),
                    Dois = LinkedIds(Cell(PublicationsColumn), DoiLink),
                    AnnounceDate = Cell(AnnounceDateColumn),
                    Keywords = Cell(KeywordsColumn),
                });
            }
            return rows;
        }

        // The publications cell is HTML (measured 2026-11-10): '<a href="https://www.ncbi.nlm.nih.gov/pubmed/42311858"
        // ...>Wu et al. (2026)</a>', '<a href="https://doi.org/10.1016/j.mcpro.2026.101658" ...>', "no publication", or
        // "Dataset with its publication pending".
        private static readonly Regex Href = new("href=\"([^\"]+)\"", RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
        private static readonly Regex PubMedLink = new(@"(?:ncbi\.nlm\.nih\.gov/pubmed/|pubmed\.ncbi\.nlm\.nih\.gov/)(\d+)", RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
        private static readonly Regex DoiLink = new(@"doi\.org/(10\.[^\s?#]+)", RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
        private static readonly Regex Tag = new("<[^>]*>", RegexOptions.CultureInvariant);

        /// <summary>The ids <paramref name="pattern"/> finds in a cell's links, in order, once each; empty for none; null for no cell.</summary>
        private static IReadOnlyList<string> LinkedIds(string cell, Regex pattern) =>
            cell == null
                ? null
                : Href.Matches(cell)
                    .Select(href => pattern.Match(Uri.UnescapeDataString(href.Groups[1].Value)))
                    .Where(m => m.Success)
                    .Select(m => m.Groups[1].Value.TrimEnd('.', ',', ';'))
                    .Distinct(StringComparer.OrdinalIgnoreCase)
                    .ToList();

        /// <summary>A cell's text with its markup removed ("Wu et al. (2026)", "no publication"), or null for no cell.</summary>
        private static string PublicationTextOf(string cell) =>
            cell == null ? null : Regex.Replace(Tag.Replace(cell, string.Empty), @"\s+", " ").Trim();

        /// <summary>Releases the <see cref="HttpClient"/> when this client created it.</summary>
        public void Dispose()
        {
            if (_disposed)
                return;
            if (_ownsHttpClient)
                _httpClient.Dispose();
            _disposed = true;
        }
    }

    /// <summary>
    /// One dataset as ProteomeCentral's search lists it. Every value is ProteomeCentral's display text, HTML entities
    /// decoded, or null when the row leaves it empty. The lab head (a person's name) is not exposed.
    /// </summary>
    public sealed class ProteomeCentralDataset
    {
        /// <summary>The ProteomeXchange accession, e.g. "PXD062174". Never null.</summary>
        public string Accession { get; init; }

        /// <summary>The title.</summary>
        public string Title { get; init; }

        /// <summary>The repository that hosts it: "PRIDE", "MassIVE", "jPOST", "iProX", "PanoramaPublic", ...</summary>
        public string Repository { get; init; }

        /// <summary>The species, as ProteomeCentral names them (several are comma-separated).</summary>
        public string Species { get; init; }

        /// <summary>ProteomeCentral's SDRF status, e.g. "curated valid SDRF", or null when it states none.</summary>
        public string Sdrf { get; init; }

        /// <summary>The "files (raw/total)" cell as ProteomeCentral gives it, or null.</summary>
        public string Files { get; init; }

        /// <summary>The instruments (several are comma-separated).</summary>
        public string Instrument { get; init; }

        /// <summary>
        /// The PubMed IDs the publications cell links to, in order. Empty when it states none ("no publication", or
        /// "Dataset with its publication pending"); null when the row has no publications cell.
        /// </summary>
        public IReadOnlyList<string> PubMedIds { get; init; }

        /// <summary>The DOIs the publications cell links to, in order; empty or null as <see cref="PubMedIds"/>.</summary>
        public IReadOnlyList<string> Dois { get; init; }

        /// <summary>
        /// The publications cell as text, its links' markup removed: "Wu et al. (2026)", "no publication", "Dataset with
        /// its publication pending". Null when the row has no publications cell.
        /// </summary>
        public string PublicationText { get; init; }

        /// <summary>The announce date, e.g. "2026-04-28".</summary>
        public string AnnounceDate { get; init; }

        /// <summary>The submitter keywords (comma-separated).</summary>
        public string Keywords { get; init; }

        /// <inheritdoc />
        public override string ToString() => $"{Accession} ({Repository}): {Title}";
    }

    /// <summary>
    /// Filters for <see cref="ProteomeCentralClient.SearchDatasetsAsync"/>. Each is a ProteomeCentral facet, measured to
    /// filter exactly to the facet's count (2026-11-10). A value must match ProteomeCentral's display name; one it does
    /// not know gives an empty result. Null leaves a filter unset.
    /// </summary>
    public sealed class ProteomeCentralSearchFilter
    {
        /// <summary>The hosting repository, e.g. "PRIDE" or "MassIVE".</summary>
        public string Repository { get; init; }

        /// <summary>
        /// The species, by ProteomeCentral's name for it, e.g. "Homo sapiens" or "Saccharomyces cerevisiae (Bakers yeast)".
        /// A taxonomy ID is refused: ProteomeCentral answers <c>species=9606</c> (or <c>NCBITaxon:9606</c>) with nothing
        /// at all, against 27,887 datasets for "Homo sapiens" (2026-11-10). MassIVE's own search takes the taxonomy ID.
        /// </summary>
        public string Species { get; init; }

        /// <summary>The instrument, e.g. "Q Exactive".</summary>
        public string Instrument { get; init; }

        /// <summary>The announce year, e.g. "2026".</summary>
        public string Year { get; init; }

        /// <summary>A submitter keyword, e.g. "DatasetType:Proteomics" (a facet; the free-text search is the keyword argument).</summary>
        public string Keywords { get; init; }

        /// <summary>The SDRF status: "No SDRF", "curated valid SDRF", "two sources SDRF", "agentic valid SDRF".</summary>
        public string Sdrf { get; init; }

        /// <summary>The file-count bucket, e.g. "0", "1-5", "6-10", "11-20".</summary>
        public string Files { get; init; }

        private static readonly Regex TaxonomyId = new(@"^(?:NCBITaxon:|taxid:)?\d+$", RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);

        /// <summary>The parameters to send, in a fixed order; a blank value, or a taxonomy ID for the species, is refused.</summary>
        internal IReadOnlyList<(string Name, string Value)> Parameters()
        {
            var parameters = new List<(string, string)>();
            void Add(string name, string value)
            {
                if (value == null)
                    return;
                if (string.IsNullOrWhiteSpace(value))
                    throw new ArgumentException($"The ProteomeCentral '{name}' filter is blank; leave it null to leave it unset.", name);
                parameters.Add((name, value));
            }

            if (Species != null && TaxonomyId.IsMatch(Species.Trim()))
                throw new ArgumentException(
                    $"ProteomeCentral filters species by name (e.g. \"Homo sapiens\"); the taxonomy ID '{Species}' would match nothing.",
                    nameof(Species));

            Add("repository", Repository);
            Add("species", Species);
            Add("instrument", Instrument);
            Add("year", Year);
            Add("keywords", Keywords);
            Add("sdrf", Sdrf);
            Add("files", Files);
            return parameters;
        }
    }
}
