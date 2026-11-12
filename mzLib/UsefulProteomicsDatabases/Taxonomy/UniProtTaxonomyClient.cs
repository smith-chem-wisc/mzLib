#nullable enable
using MzLibUtil;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Net.Http;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;

namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>A taxonomy built from UniProt: the tree, and the proteomes that were read to build it.</summary>
    /// <param name="Tree">The organisms of the proteomes read, and every ancestor of each.</param>
    /// <param name="Proteomes">
    /// The proteomes read, each with the taxon it belongs to. Every taxon named is in <paramref name="Tree"/>, except
    /// those of <paramref name="ProteomesWithoutLineage"/>.
    /// </param>
    /// <param name="ProteomesWithoutLineage">
    /// The identifiers, in ordinal order, of proteomes whose taxon UniProt names by number only, with no name and no
    /// lineage. These are taxa NCBI has since merged or withdrawn. Such a proteome is in <paramref name="Proteomes"/>
    /// with its taxon identifier, but that taxon is not in <paramref name="Tree"/>, so asking the tree about it
    /// answers <see cref="TaxonRankStatus.TaxonNotInTree"/> rather than inventing a place for it.
    /// </param>
    /// <param name="ProteomesNotFound">
    /// The identifiers, in ordinal order, that <see cref="UniProtTaxonomyClient.FetchProteomesAsync"/> was asked for
    /// and UniProt did not return. Always empty for a fetch by query.
    /// </param>
    public sealed record UniProtTaxonomySnapshot(TaxonomyTree Tree, ProteomeTaxonMap Proteomes,
        IReadOnlyList<string> ProteomesWithoutLineage, IReadOnlyList<string> ProteomesNotFound);

    /// <summary>
    /// Builds a <see cref="TaxonomyTree"/> and a <see cref="ProteomeTaxonMap"/> from UniProt's proteome catalogue.
    /// </summary>
    /// <remarks>
    /// <para>
    /// UniProt describes each proteome with its organism's taxon and that taxon's whole lineage, every ancestor with
    /// its identifier, name and rank. Reading the proteomes that match a query therefore yields both which taxon each
    /// proteome belongs to and the tree above all of them. Nothing is inferred from names.
    /// </para>
    /// <para>
    /// One thing the proteome listing leaves out is the rank of the organism itself: it gives Homo sapiens' ancestors
    /// their ranks but not Homo sapiens its own. Those are read in a second step from UniProt's taxonomy endpoint, a
    /// few hundred organisms per request, so that "the species of this proteome" has an answer when the proteome's
    /// organism is the species.
    /// </para>
    /// <para>
    /// The result is as large as the query: the reference proteomes (<see cref="ReferenceProteomesQuery"/>, the
    /// default, tens of thousands) take a couple of minutes; every proteome UniProt holds is over a million. mzLib
    /// ships the reference set already built, as <see cref="TaxonomyTree.LoadPinned"/>; this client is how that
    /// snapshot is refreshed and how a caller builds a different one.
    /// </para>
    /// <para>
    /// Failures follow the split the other repository clients use: an <see cref="MzLibException"/> when UniProt
    /// answered with something this client cannot use, an <see cref="HttpRequestException"/> when it could not be
    /// reached or failed. Each page is retried on a transient failure before the fetch gives up.
    /// </para>
    /// </remarks>
    public sealed class UniProtTaxonomyClient : IDisposable
    {
        /// <summary>The query selecting UniProt's reference proteomes.</summary>
        public const string ReferenceProteomesQuery = "reference:true";

        /// <summary>The endpoint the proteomes are read from, recorded as the snapshot's source.</summary>
        public const string ProteomeSearchUrl = "https://rest.uniprot.org/proteomes/search";

        /// <summary>The endpoint the organisms' own ranks are read from.</summary>
        public const string TaxonomyUrl = "https://rest.uniprot.org/taxonomy/taxonIds";

        /// <summary>The most proteomes UniProt serves in one page.</summary>
        public const int MaxPageSize = 500;

        // Identifiers per taxonomy request. They travel in the URL, so the batch is bounded by its length:
        // 200 seven-digit identifiers is about 1,600 characters.
        private const int TaxaPerRequest = 200;

        // Proteome identifiers per query when fetching a named set; "upid:UP000005640 OR " is 20 characters.
        private const int ProteomesPerQuery = 100;

        private static readonly Regex ProteomeIdPattern = new(@"^UP\d{9,}$", RegexOptions.Compiled);

        // A fetch that has not finished after this many pages is not paging: UniProt's whole catalogue is under
        // 3,000 pages of 500.
        private const int MaxPages = 20000;

        private static readonly Regex NextLink = new(@"<([^>]+)>\s*;\s*rel=""next""", RegexOptions.Compiled);

        private readonly HttpClient _httpClient;
        private readonly bool _ownsHttpClient;
        private readonly DownloadEngine _engine;

        /// <summary>Creates a client with its own connection to UniProt.</summary>
        public UniProtTaxonomyClient()
            : this(new HttpClient { Timeout = TimeSpan.FromSeconds(100) }, ownsHttpClient: true)
        {
        }

        /// <summary>Creates a client over a caller's <see cref="HttpClient"/>, which the client does not dispose.</summary>
        public UniProtTaxonomyClient(HttpClient httpClient)
            : this(httpClient, ownsHttpClient: false)
        {
        }

        /// <summary>Internal so tests can check that an owned connection is disposed and a borrowed one is not.</summary>
        internal UniProtTaxonomyClient(HttpClient httpClient, bool ownsHttpClient)
        {
            _httpClient = httpClient ?? throw new ArgumentNullException(nameof(httpClient));
            _ownsHttpClient = ownsHttpClient;
            _engine = new DownloadEngine(_httpClient, "UniProt");
        }

        /// <summary>The engine behind every request; tests shorten its retry delays.</summary>
        internal DownloadEngine Engine => _engine;

        /// <summary>
        /// Reads every proteome matching <paramref name="query"/> and builds the taxonomy above them.
        /// </summary>
        /// <param name="query">A UniProt proteome query, e.g. <see cref="ReferenceProteomesQuery"/> or
        /// "taxonomy_id:10239" for every viral proteome.</param>
        /// <param name="pageSize">Proteomes per request, from 1 to <see cref="MaxPageSize"/>.</param>
        /// <param name="cancellationToken">Cancels the fetch between or during pages.</param>
        /// <exception cref="ArgumentException"><paramref name="query"/> is blank.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="pageSize"/> is out of range.</exception>
        /// <exception cref="MzLibException">UniProt answered with something that is not a proteome listing, matched no
        /// proteome, described a proteome without an identifier or taxon, or gave one taxon two different parents.</exception>
        /// <exception cref="HttpRequestException">UniProt could not be reached, or failed, after retries.</exception>
        public async Task<UniProtTaxonomySnapshot> FetchAsync(string query = ReferenceProteomesQuery, int pageSize = MaxPageSize,
            CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(query))
            {
                throw new ArgumentException("A proteome query is required.", nameof(query));
            }
            if (pageSize < 1 || pageSize > MaxPageSize)
            {
                throw new ArgumentOutOfRangeException(nameof(pageSize), pageSize, $"The page size must be from 1 to {MaxPageSize}.");
            }

            var read = new Listing();
            await ReadQueryAsync(query, pageSize, read, cancellationToken).ConfigureAwait(false);
            if (read.Proteomes.Count == 0)
            {
                throw new MzLibException($"UniProt has no proteome matching '{query}'.");
            }
            return await BuildAsync(read, query, Array.Empty<string>(), cancellationToken).ConfigureAwait(false);
        }

        /// <summary>
        /// Reads a named set of proteomes and builds the taxonomy above them: the tree for exactly the databases a
        /// search uses, whether or not they are reference proteomes.
        /// </summary>
        /// <remarks>
        /// The tree mzLib ships covers reference proteomes only. A folder of downloaded proteomes usually holds a few
        /// that are not, or that have lost reference status since; this places all of them, from one release, without
        /// overlaying one release's taxa on another's. An identifier UniProt does not return is listed in
        /// <see cref="UniProtTaxonomySnapshot.ProteomesNotFound"/> rather than failing the fetch.
        /// </remarks>
        /// <param name="proteomeIds">UniProt proteome identifiers, e.g. "UP000005640". Repeats are read once.</param>
        /// <param name="cancellationToken">Cancels the fetch between or during requests.</param>
        /// <exception cref="ArgumentException">There are no identifiers, or one is not of the form "UP" followed by digits.</exception>
        /// <exception cref="MzLibException">UniProt returned none of the proteomes, or answered with something that is
        /// not a proteome listing.</exception>
        /// <exception cref="HttpRequestException">UniProt could not be reached, or failed, after retries.</exception>
        public async Task<UniProtTaxonomySnapshot> FetchProteomesAsync(IEnumerable<string> proteomeIds,
            CancellationToken cancellationToken = default)
        {
            ArgumentNullException.ThrowIfNull(proteomeIds);
            List<string> ids = proteomeIds.Distinct(StringComparer.Ordinal).OrderBy(id => id, StringComparer.Ordinal).ToList();
            if (ids.Count == 0)
            {
                throw new ArgumentException("At least one proteome identifier is required.", nameof(proteomeIds));
            }
            foreach (string id in ids)
            {
                if (id == null || !ProteomeIdPattern.IsMatch(id))
                {
                    throw new ArgumentException($"'{id}' is not a UniProt proteome identifier such as UP000005640.", nameof(proteomeIds));
                }
            }

            var read = new Listing();
            for (int start = 0; start < ids.Count; start += ProteomesPerQuery)
            {
                string query = string.Join(" OR ", ids.Skip(start).Take(ProteomesPerQuery).Select(id => "upid:" + id));
                await ReadQueryAsync(query, MaxPageSize, read, cancellationToken).ConfigureAwait(false);
            }

            List<string> notFound = ids.Where(id => !read.Proteomes.ContainsKey(id)).ToList();
            if (read.Proteomes.Count == 0)
            {
                throw new MzLibException($"UniProt returned none of the {ids.Count} proteomes asked for.");
            }
            string described = ids.Count == 1 ? "upid:" + ids[0] : $"{ids.Count} proteomes by identifier";
            return await BuildAsync(read, described, notFound, cancellationToken).ConfigureAwait(false);
        }

        /// <summary>What has been read so far, across the pages of one or more queries.</summary>
        private sealed class Listing
        {
            public Dictionary<long, TaxonomyNode> Nodes { get; } = new();
            public Dictionary<string, ProteomeTaxon> Proteomes { get; } = new(StringComparer.Ordinal);
            public List<string> WithoutLineage { get; } = new();
            public HashSet<string> Requested { get; } = new(StringComparer.Ordinal);
            public string? Release { get; set; }
        }

        /// <summary>Reads every page of one query into <paramref name="read"/>.</summary>
        private async Task ReadQueryAsync(string query, int pageSize, Listing read, CancellationToken cancellationToken)
        {
            string? requestUri = ProteomeSearchUrl + "?query=" + Uri.EscapeDataString(query)
                + "&size=" + pageSize.ToString(CultureInfo.InvariantCulture);

            for (int page = 0; requestUri != null; page++)
            {
                cancellationToken.ThrowIfCancellationRequested();
                if (page >= MaxPages || !read.Requested.Add(requestUri))
                {
                    throw new MzLibException(
                        $"UniProt's proteome listing for '{query}' did not end: page {page + 1} repeats an earlier request or exceeds {MaxPages} pages.");
                }

                string thisRequest = requestUri;
                using HttpResponseMessage response = await _engine.GetSuccessAsync(thisRequest,
                    r => $"UniProt proteome request failed with status {(int)r.StatusCode} {r.ReasonPhrase} for query '{query}'.",
                    cancellationToken).ConfigureAwait(false);

                read.Release ??= HeaderValue(response, "x-uniprot-release");
                string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

                JArray results = ParseResults(content, query);
                foreach (JToken result in results)
                {
                    ReadProteome(result, read.Nodes, read.Proteomes, read.WithoutLineage);
                }

                // An empty page ends the listing even if it carries a next link, so a server that keeps
                // offering one cannot hold the fetch open.
                requestUri = results.Count == 0 ? null : NextLinkOf(response);
            }
        }

        /// <summary>Fills in the organisms' ranks and turns what was read into a tree and a map.</summary>
        private async Task<UniProtTaxonomySnapshot> BuildAsync(Listing read, string query, IReadOnlyList<string> notFound,
            CancellationToken cancellationToken)
        {
            await FillOrganismRanksAsync(read.Nodes, cancellationToken).ConfigureAwait(false);

            var provenance = new TaxonomyProvenance(ProteomeSearchUrl, read.Release, query);
            TaxonomyTree tree;
            try
            {
                tree = TaxonomyTree.FromNodes(read.Nodes.Values, provenance);
            }
            catch (System.IO.InvalidDataException e)
            {
                throw new MzLibException($"The lineages UniProt served for '{query}' do not form a tree: {e.Message}", e);
            }
            read.WithoutLineage.Sort(StringComparer.Ordinal);
            return new UniProtTaxonomySnapshot(tree, ProteomeTaxonMap.FromProteomes(read.Proteomes.Values, provenance),
                read.WithoutLineage, notFound);
        }

        /// <summary>
        /// Gives each taxon that still has no rank the one UniProt's taxonomy endpoint states for it. These are the
        /// proteomes' own organisms, which the proteome listing names without a rank. A taxon the endpoint does not
        /// return keeps its empty rank: the source stated none.
        /// </summary>
        private async Task FillOrganismRanksAsync(Dictionary<long, TaxonomyNode> nodes, CancellationToken cancellationToken)
        {
            List<long> unranked = nodes.Values.Where(n => n.Rank.Length == 0).Select(n => n.TaxonId).OrderBy(id => id).ToList();
            for (int start = 0; start < unranked.Count; start += TaxaPerRequest)
            {
                cancellationToken.ThrowIfCancellationRequested();
                string ids = string.Join(",", unranked.Skip(start).Take(TaxaPerRequest)
                    .Select(id => id.ToString(CultureInfo.InvariantCulture)));
                string requestUri = TaxonomyUrl + "/" + ids + "?fields=id,rank,parent&size=" + MaxPageSize.ToString(CultureInfo.InvariantCulture);

                using HttpResponseMessage response = await _engine.GetSuccessAsync(requestUri,
                    r => $"UniProt taxonomy request failed with status {(int)r.StatusCode} {r.ReasonPhrase}.",
                    cancellationToken).ConfigureAwait(false);
                string content = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);

                foreach (JToken taxon in ParseResults(content, "taxonomy/taxonIds"))
                {
                    if (Whole(taxon, "taxonId") is not long taxonId || !nodes.TryGetValue(taxonId, out var node))
                    {
                        continue;
                    }
                    long? statedParent = taxon["parent"] is JObject parent ? Whole(parent, "taxonId") : null;
                    if (statedParent != null && node.ParentTaxonId != null && statedParent != node.ParentTaxonId)
                    {
                        throw new MzLibException(
                            $"UniProt's taxonomy places taxon {taxonId} ({node.ScientificName}) under {Describe(statedParent)}, " +
                            $"but its proteome's lineage ends at {Describe(node.ParentTaxonId)}.");
                    }
                    nodes[taxonId] = node with { Rank = Text(taxon, "rank") ?? "" };
                }
            }
        }

        /// <inheritdoc/>
        public void Dispose()
        {
            if (_ownsHttpClient)
            {
                _httpClient.Dispose();
            }
        }

        private static JArray ParseResults(string content, string query)
        {
            JToken? results;
            try
            {
                results = JObject.Parse(content)["results"];
            }
            catch (JsonException e)
            {
                throw new MzLibException($"UniProt's answer to '{query}' is not JSON: {e.Message}", e);
            }
            return results as JArray
                ?? throw new MzLibException($"UniProt's answer to '{query}' has no 'results' list.");
        }

        /// <summary>
        /// Adds one proteome, its organism and its organism's ancestors. The lineage is listed root first, so each
        /// entry's parent is the entry before it, and the organism's parent is the last entry.
        /// </summary>
        private static void ReadProteome(JToken result, Dictionary<long, TaxonomyNode> nodes, Dictionary<string, ProteomeTaxon> proteomes,
            List<string> withoutLineage)
        {
            string? proteomeId = Text(result, "id");
            if (string.IsNullOrEmpty(proteomeId))
            {
                throw new MzLibException("UniProt listed a proteome with no identifier.");
            }
            if (result["taxonomy"] is not JObject organism || Whole(organism, "taxonId") is not long organismId)
            {
                throw new MzLibException($"UniProt listed proteome {proteomeId} with no taxon.");
            }

            long? parent = null;
            bool organismInLineage = false;
            if (result["taxonLineage"] is JArray lineage)
            {
                foreach (JToken entry in lineage)
                {
                    if (Whole(entry, "taxonId") is not long taxonId)
                    {
                        throw new MzLibException($"UniProt listed proteome {proteomeId} with a lineage entry that has no taxon identifier.");
                    }
                    Add(nodes, NodeOf(entry, taxonId, parent), proteomeId);
                    organismInLineage |= taxonId == organismId;
                    parent = taxonId;
                }
            }
            // A taxon NCBI has merged or withdrawn comes with its number and nothing else. There is no name to
            // show and no parent to hang it from, so it gets no node; the proteome is still recorded below.
            bool described = parent != null || !string.IsNullOrEmpty(Text(organism, "scientificName"));
            if (!described)
            {
                withoutLineage.Add(proteomeId);
            }
            else if (!organismInLineage)
            {
                Add(nodes, NodeOf(organism, organismId, parent), proteomeId);
            }

            var proteome = new ProteomeTaxon(proteomeId, organismId,
                Text(organism, "scientificName") ?? "",
                Whole(result, "proteinCount"),
                Text(result, "proteomeType"));
            if (!proteomes.TryAdd(proteomeId, proteome))
            {
                throw new MzLibException($"UniProt listed proteome {proteomeId} twice.");
            }
        }

        private static TaxonomyNode NodeOf(JToken source, long taxonId, long? parent) =>
            new(taxonId, parent,
                Text(source, "rank") ?? "",
                Text(source, "scientificName") ?? "",
                Text(source, "commonName"),
                source["hidden"]?.Type == JTokenType.Boolean && source.Value<bool>("hidden"));

        /// <summary>
        /// Keeps the first description of a taxon. A later one must put it under the same parent: two proteomes
        /// that disagree about where a taxon sits cannot both be placed in one tree.
        /// </summary>
        private static void Add(Dictionary<long, TaxonomyNode> nodes, TaxonomyNode node, string proteomeId)
        {
            if (!nodes.TryGetValue(node.TaxonId, out var existing))
            {
                nodes[node.TaxonId] = node;
                return;
            }
            if (existing.ParentTaxonId != node.ParentTaxonId)
            {
                throw new MzLibException(
                    $"UniProt places taxon {node.TaxonId} ({node.ScientificName}) under {Describe(node.ParentTaxonId)} in the lineage of " +
                    $"proteome {proteomeId}, but under {Describe(existing.ParentTaxonId)} in an earlier proteome's.");
            }
            // The organism of a proteome is named without a rank; a lineage that passes through it states one.
            if (existing.Rank.Length == 0 && node.Rank.Length != 0)
            {
                nodes[node.TaxonId] = existing with { Rank = node.Rank, IsHidden = node.IsHidden };
            }
        }

        private static string Describe(long? taxonId) =>
            taxonId is long id ? id.ToString(CultureInfo.InvariantCulture) : "no parent";

        private static string? Text(JToken source, string property)
        {
            JToken? value = source[property];
            return value == null || value.Type == JTokenType.Null ? null : value.ToString();
        }

        private static long? Whole(JToken source, string property)
        {
            JToken? value = source[property];
            return value != null && value.Type == JTokenType.Integer ? value.Value<long>() : null;
        }

        private static string? HeaderValue(HttpResponseMessage response, string name) =>
            response.Headers.TryGetValues(name, out var values) ? values.FirstOrDefault() : null;

        /// <summary>The URL of the next page, from the response's "Link: &lt;url&gt;; rel="next"" header; null on the last page.</summary>
        private static string? NextLinkOf(HttpResponseMessage response)
        {
            if (!response.Headers.TryGetValues("Link", out var links))
            {
                return null;
            }
            foreach (string link in links)
            {
                Match match = NextLink.Match(link);
                if (match.Success)
                {
                    return match.Groups[1].Value;
                }
            }
            return null;
        }
    }
}
