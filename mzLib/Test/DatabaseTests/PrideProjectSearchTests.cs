using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideProjectSearchTests
{
    // ---- test doubles -------------------------------------------------------

    /// <summary>An HttpMessageHandler that returns a caller-supplied response and records request URIs.</summary>
    private sealed class StubHandler : HttpMessageHandler
    {
        private readonly Func<HttpRequestMessage, HttpResponseMessage> _responder;
        public List<string> RequestedUris { get; } = new();

        public StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) => _responder = responder;

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            RequestedUris.Add(request.RequestUri.ToString());
            return Task.FromResult(_responder(request));
        }
    }

    private static HttpResponseMessage JsonResponse(string body, HttpStatusCode status = HttpStatusCode.OK, long? totalRecords = null)
    {
        var response = new HttpResponseMessage(status) { Content = new StringContent(body, Encoding.UTF8, "application/json") };
        if (totalRecords.HasValue)
            response.Headers.TryAddWithoutValidation("total_records", totalRecords.Value.ToString());
        return response;
    }

    /// <summary>
    /// One PRIDE search hit in the real v3 wire shape, captured live from
    /// <c>search/projects?keyword=hepatocyte</c> on 2026-08-21 (PXD068742).
    /// </summary>
    /// <remarks>
    /// TRIMMED, NOT EDITED: collections are cut to their first two entries and free text is truncated
    /// with an ellipsis. No field was added, removed, renamed or given a different value. That
    /// distinction matters — an earlier revision of this fixture carried an <c>sdrf</c> member and
    /// invented <c>yearlyDownloads</c> counts, neither of which the server sends for this project, so
    /// two tests pinned numbers PRIDE does not serve. Re-capture rather than hand-edit.
    /// <para>
    /// Note what the shape shows: every collection is an array of plain STRINGS — no CvParam, no
    /// contact objects, no structured references — which is the whole reason this endpoint needs its
    /// own DTO. Also the bare calendar dates, the "organismsPart" spelling (the metadata endpoint says
    /// "organismParts"), the dynamic "highlights" map carrying &lt;em&gt; markup, and the pre-mangled
    /// reference string with its literal "null" segments.
    /// </para>
    /// <para>
    /// This project reports NO <c>sdrf</c> member at all, and its <c>doi</c>, <c>projectTags</c> and
    /// <c>diseases</c> are empty. Those are exercised against
    /// <see cref="PopulatedOptionalsHitJson"/> instead — asserting them empty here would pass whether
    /// the property bound correctly or bound nothing at all.
    /// </para>
    /// </remarks>
    private const string SearchHitJson =
        """
        {
          "accession": "PXD068742",
          "title": "Culturing of primary mouse hepatocytes is associated with major change...",
          "projectDescription": "Primary hepatocytes are a commonly used in vitro model for studying li...",
          "dataProcessingProtocol": "The FragPipe-generated master spectral library was used in concert wit...",
          "sampleProcessingProtocol": "Samples were homogenized in 4% SDS buffer (0.1 M TrisHCl and 4% SDS),...",
          "projectTags": [],
          "keywords": [
            "Cell culture",
            "Hepatocytes"
          ],
          "doi": "",
          "otherOmicsLinks": [
            "geo:GSE280301",
            "geo:GSE173406"
          ],
          "submissionType": "PARTIAL",
          "publicationDate": "2025-09-25",
          "updatedDate": "2025-09-23",
          "submissionDate": "2025-09-23",
          "downloadCount": 111,
          "avgDownloadsPerFile": 3.9642856,
          "botCount": 52,
          "hubCount": 47,
          "organicCount": 12,
          "percentile": 5,
          "yearlyDownloads": [
            {
              "year": "2025",
              "count": 33
            },
            {
              "year": "2026",
              "count": 78
            }
          ],
          "submitters": [
            "Ben Stocks"
          ],
          "labPIs": [
            "Jonas T. Treebak"
          ],
          "affiliations": [
            "Novo Nordisk Foundation Center for Basic Metabolic Research"
          ],
          "instruments": [
            "timsTOF Pro 2"
          ],
          "softwares": [
            "FragPipe",
            "DIA-NN"
          ],
          "quantificationMethods": [
            "diaPASEF"
          ],
          "sampleAttributes": [
            "hepatocyte",
            "liver"
          ],
          "organisms": [
            "Mus musculus (mouse)"
          ],
          "organismsPart": [
            "Hepatocyte",
            "Liver"
          ],
          "diseases": [],
          "references": [
            "null--pubMed:0--doi: 10.1097/HC9.0000000000000795"
          ],
          "experimentTypes": [
            "Data-independent acquisition"
          ],
          "projectFileNames": [
            "cs_quant.pg_matrix.tsv",
            "20210930_PRI_TIMS4_PRIEvo1_P094_R0275_44min_S12_S2-D12_1_2797.d.rar"
          ],
          "highlights": {
            "projectDescription": [
              "Cultures contained ~10% non-<em>hepatocytes</em>, including stellate- and dedifferentiat..."
            ]
          }
        }
        """;

    /// <summary>
    /// A second live hit (PXD080181, same capture date) carrying optional members the primary fixture
    /// leaves empty: a populated <c>sdrf</c>, <c>projectTags</c> and <c>diseases</c>. Trimmed on the
    /// same terms — entries dropped and prose truncated, nothing invented.
    /// </summary>
    /// <remarks>
    /// It exists so those members are pinned by an assertion that can actually FAIL. Note <c>sdrf</c>
    /// is not a file or a filename: it is a flattened bag of the SDRF's term values, which is why the
    /// DTO documents it that way.
    /// </remarks>
    private const string PopulatedOptionalsHitJson =
        """
        {
          "accession": "PXD080181",
          "title": "Platycodi Radix Extract Prevents Hepatic Steatosis by Enhancing Bile A...",
          "projectDescription": "This project investigated the molecular mechanism by which *Platycodi...",
          "dataProcessingProtocol": "The original MASCOT 2.4 result files were no longer available, so the...",
          "sampleProcessingProtocol": "Five-week-old wild-type male C57BL/6J mice (Orient Bio, Seongnam, Kore...",
          "projectTags": [
            "Liver (b/d-hpp)",
            "Human proteome project"
          ],
          "keywords": [
            "Obesity",
            "Nafld"
          ],
          "doi": "",
          "submissionType": "PARTIAL",
          "publicationDate": "2026-07-04",
          "updatedDate": "2026-06-25",
          "submissionDate": "2026-06-25",
          "downloadCount": 0,
          "avgDownloadsPerFile": 0.0,
          "percentile": 0,
          "submitters": [
            "Wooyoung Kim"
          ],
          "labPIs": [
            "Edmond Changkyun Park"
          ],
          "affiliations": [
            "Korean Basic Science Instittue"
          ],
          "instruments": [
            "LTQ Orbitrap Velos"
          ],
          "softwares": [],
          "quantificationMethods": [],
          "sampleAttributes": [
            "liver",
            "Mus musculus (Mouse)"
          ],
          "organisms": [
            "Mus musculus",
            "Mus musculus (mouse)"
          ],
          "organismsPart": [
            "Liver"
          ],
          "diseases": [
            "Normal",
            "Non-alcoholic fatty liver disease"
          ],
          "references": [
            "Kim W, Baek WH, Yun SH, Lee H, Kim MJ, Lee SY, Kim GH, Kim SI, Jeong H..."
          ],
          "experimentTypes": [
            "Data-dependent acquisition"
          ],
          "sdrf": "Nt=label free sample;ac=ms:1002038 Nt=trypsin;ac=ms:1001251 Mus muscul...",
          "projectFileNames": [
            "150609_YSH_LIVER_HFD_10_06.raw",
            "150811_YSH_LIVER_NCD_10_07.raw"
          ],
          "highlights": {
            "sdrf": [
              "...ac=ms:1001742 Non-alcoholic fatty liver <em>disease</em>"
            ]
          }
        }
        """;

    /// <summary>A minimal hit carrying only the accession — enough to identify it while paging.</summary>
    private static string HitJson(string accession) => $$"""{ "accession": "{{accession}}" }""";

    private static string Array(params string[] objects) => "[" + string.Join(",", objects) + "]";

    // ---- deserialization ----------------------------------------------------

    [Test]
    public async Task SearchProjectsAsync_RealWireShape_PopulatesEveryField()
    {
        var handler = new StubHandler(_ => JsonResponse(Array(SearchHitJson)));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("liver");

        Assert.That(hits.Count, Is.EqualTo(1));
        PrideProjectSearchResult hit = hits[0];
        Assert.Multiple(() =>
        {
            Assert.That(hit.Accession, Is.EqualTo("PXD068742"));
            Assert.That(hit.Title, Does.StartWith("Culturing of primary mouse hepatocytes"));
            Assert.That(hit.ProjectDescription, Does.StartWith("Primary hepatocytes"));
            Assert.That(hit.DataProcessingProtocol, Does.Contain("FragPipe"));
            Assert.That(hit.SampleProcessingProtocol, Does.Contain("SDS buffer"));
            Assert.That(hit.SubmissionType, Is.EqualTo("PARTIAL"));
            // Doi, Sdrf, ProjectTags and Diseases are deliberately NOT asserted here. They are empty
            // (Sdrf absent outright) for this project, so "Is.Empty" would hold whether the property
            // bound correctly or bound nothing at all — the exact silent-bind failure OrganismParts
            // has its own test for. They are pinned against PopulatedOptionalsHitJson instead.

            // Bare calendar dates: no time component, and no offset attached from the local machine.
            Assert.That(hit.SubmissionDate, Is.EqualTo(new DateTime(2025, 9, 23)));
            Assert.That(hit.PublicationDate, Is.EqualTo(new DateTime(2025, 9, 25)));
            Assert.That(hit.UpdatedDate, Is.EqualTo(new DateTime(2025, 9, 23)));
            Assert.That(hit.SubmissionDate.TimeOfDay, Is.EqualTo(TimeSpan.Zero));

            // Every controlled-vocabulary collection arrives FLATTENED to display strings. This is the
            // assertion that would fail first if PRIDE ever served CvParam objects here, which is the
            // drift that would invalidate the whole DTO.
            Assert.That(hit.Instruments, Is.EqualTo(new[] { "timsTOF Pro 2" }));
            Assert.That(hit.Softwares, Is.EqualTo(new[] { "FragPipe", "DIA-NN" }));
            Assert.That(hit.QuantificationMethods, Is.EqualTo(new[] { "diaPASEF" }));
            Assert.That(hit.Organisms, Is.EqualTo(new[] { "Mus musculus (mouse)" }));
            Assert.That(hit.ExperimentTypes, Is.EqualTo(new[] { "Data-independent acquisition" }));
            Assert.That(hit.SampleAttributes, Is.EqualTo(new[] { "hepatocyte", "liver" }));
            Assert.That(hit.Keywords, Is.EqualTo(new[] { "Cell culture", "Hepatocytes" }));
            Assert.That(hit.Submitters, Is.EqualTo(new[] { "Ben Stocks" }));
            Assert.That(hit.LabPIs, Is.EqualTo(new[] { "Jonas T. Treebak" }));
            Assert.That(hit.Affiliations, Is.EqualTo(new[] { "Novo Nordisk Foundation Center for Basic Metabolic Research" }));
            Assert.That(hit.References, Is.EqualTo(new[] { "null--pubMed:0--doi: 10.1097/HC9.0000000000000795" }));
            Assert.That(hit.ProjectFileNames, Is.EqualTo(new[] { "cs_quant.pg_matrix.tsv", "20210930_PRI_TIMS4_PRIEvo1_P094_R0275_44min_S12_S2-D12_1_2797.d.rar" }));
            Assert.That(hit.OtherOmicsLinks, Is.EqualTo(new[] { "geo:GSE280301", "geo:GSE173406" }));

            Assert.That(hit.DownloadCount, Is.EqualTo(111));
            Assert.That(hit.AvgDownloadsPerFile, Is.EqualTo(3.9642856).Within(1e-7));
            Assert.That(hit.Percentile, Is.EqualTo(5));
            Assert.That(hit.BotCount, Is.EqualTo(52));
            Assert.That(hit.HubCount, Is.EqualTo(47));
            Assert.That(hit.OrganicCount, Is.EqualTo(12));

            Assert.That(hit.YearlyDownloads.Select(y => y.Year), Is.EqualTo(new[] { "2025", "2026" }));
            Assert.That(hit.YearlyDownloads.Select(y => y.Count), Is.EqualTo(new long[] { 33, 78 }));

            // Membership, not order: PRIDE serialises `highlights` from an unordered map, so its key
            // order varies between identical requests (observed live 2026-08-21). Safe against this
            // single-key stub, but an order-sensitive assertion here would invite a flaky live one.
            Assert.That(hit.Highlights.Keys, Is.EquivalentTo(new[] { "projectDescription" }));
            Assert.That(hit.Highlights["projectDescription"].Single(), Does.Contain("<em>hepatocytes</em>"));
        });
    }

    /// <summary>
    /// The optional members the primary fixture leaves empty, asserted against a hit that actually
    /// carries them. Without this they would be pinned only by "Is.Empty" assertions that pass
    /// whether the property binds or not — the same silent-bind failure OrganismParts guards against.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_OptionalMembers_BindWhenTheServerSendsThem()
    {
        var handler = new StubHandler(_ => JsonResponse(Array(PopulatedOptionalsHitJson)));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideProjectSearchResult hit = (await client.SearchProjectsAsync("liver")).Single();

        Assert.Multiple(() =>
        {
            Assert.That(hit.Accession, Is.EqualTo("PXD080181"));
            Assert.That(hit.ProjectTags, Is.EqualTo(new[] { "Liver (b/d-hpp)", "Human proteome project" }));
            Assert.That(hit.Diseases, Is.EqualTo(new[] { "Normal", "Non-alcoholic fatty liver disease" }));
            // Not a filename and not a URL — a flattened bag of the SDRF's term values.
            Assert.That(hit.Sdrf, Does.StartWith("Nt=label free sample;ac=ms:1002038"));
        });
    }

    /// <summary>
    /// Doi gets its own minimal probe because neither live fixture carries one — only 77 of 421
    /// sampled hits had a DOI, and both captured projects fall in the majority that do not.
    /// </summary>
    /// <remarks>
    /// This JSON is SYNTHETIC and says so, unlike the two captured fixtures above. It asserts one
    /// thing only: that the property binds the wire name. That is worth pinning because `doi` is one
    /// of ten members that bind through Newtonsoft's case-insensitive fallback rather than an exact
    /// match, so nothing but a test would notice if that fallback ever went away.
    /// </remarks>
    [Test]
    public async Task SearchProjectsAsync_Doi_BindsFromTheWire()
    {
        const string hitWithDoi = """{ "accession": "PXD000001", "doi": "10.6019/PXD000001" }""";
        var handler = new StubHandler(_ => JsonResponse(Array(hitWithDoi)));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideProjectSearchResult hit = (await client.SearchProjectsAsync("liver")).Single();

        Assert.That(hit.Doi, Is.EqualTo("10.6019/PXD000001"));
    }

    // ---- mid-fetch drift ----------------------------------------------------

    /// <summary>
    /// PRIDE pages a LIVE index, so a result set can change between two requests of one fetch. When a
    /// project is published mid-fetch every later record shifts down a position, and a record that was
    /// on page 1 is served again on page 2 — verified live on 2026-08-21, where a project leaving the
    /// "liver" set slid the whole window by one.
    /// </summary>
    /// <remarks>
    /// The shared pager cannot defend this: it refuses to treat any FIELD as a record's identity,
    /// because the file manifest it was written for has no unique key. A search hit does — accession
    /// is exactly that — so the reasoning does not carry over, and search dedups on it. The mirror
    /// case (a record SKIPPED because the window slid the other way) is not fixable without a stable
    /// cursor the API does not offer, and is documented on the method rather than claimed away.
    /// </remarks>
    [Test]
    public async Task SearchProjectsAsync_ResultSetGrowsMidFetch_DoesNotReturnADuplicateProject()
    {
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            // page 0 of the set as it was; page 1 after a new project pushed everything down one, so
            // PXD000002 is served a second time.
            if (uri.Contains("page=0")) return JsonResponse(Array(HitJson("PXD000001"), HitJson("PXD000002")), totalRecords: 4);
            if (uri.Contains("page=1")) return JsonResponse(Array(HitJson("PXD000002"), HitJson("PXD000003")), totalRecords: 4);
            return JsonResponse("[ ]", totalRecords: 4);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("liver", pageSize: 2);

        Assert.That(hits.Select(h => h.Accession),
            Is.EqualTo(new[] { "PXD000001", "PXD000002", "PXD000003" }),
            "a project served on two pages must be returned once");
    }

    /// <summary>
    /// The dedup must not mistake DISTINCT projects for repeats — it keys on accession precisely
    /// because that is a real identity, so nothing else about a hit may affect whether it survives.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_DistinctProjects_AreAllKept()
    {
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("page=0")) return JsonResponse(Array(HitJson("PXD000001"), HitJson("PXD000002")), totalRecords: 4);
            if (uri.Contains("page=1")) return JsonResponse(Array(HitJson("PXD000003"), HitJson("PXD000004")), totalRecords: 4);
            return JsonResponse("[ ]", totalRecords: 4);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("liver", pageSize: 2);

        Assert.That(hits.Select(h => h.Accession),
            Is.EqualTo(new[] { "PXD000001", "PXD000002", "PXD000003", "PXD000004" }));
    }

    /// <summary>
    /// The one member that cannot bind by convention. PRIDE spells the field "organismsPart" here and
    /// "organismParts" on projects/{accession}; the DTO keeps the metadata spelling so the two agree,
    /// which only works because of the explicit JsonProperty. Without it the property binds nothing
    /// and silently returns an empty list — a wrong answer, not an error — so it gets its own test.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_OrganismsPartWireName_BindsToOrganismParts()
    {
        var handler = new StubHandler(_ => JsonResponse(Array(SearchHitJson)));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("liver");

        Assert.That(hits[0].OrganismParts, Is.EqualTo(new[] { "Hepatocyte", "Liver" }));
    }

    /// <summary>
    /// Members PRIDE omits arrive as empty collections and zeros, never null — the guarantee the DTO's
    /// remarks make. Several search fields are genuinely sparse (yearlyDownloads was absent from 61 of
    /// 100 hits in the capture sample), so a caller dereferencing one is the ordinary case.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_OmittedMembers_AreEmptyNotNull()
    {
        var handler = new StubHandler(_ => JsonResponse(Array(HitJson("PXD000001"))));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideProjectSearchResult hit = (await client.SearchProjectsAsync("liver")).Single();

        Assert.Multiple(() =>
        {
            Assert.That(hit.Accession, Is.EqualTo("PXD000001"));
            Assert.That(hit.Title, Is.Empty);
            Assert.That(hit.Keywords, Is.Empty);
            Assert.That(hit.Instruments, Is.Empty);
            Assert.That(hit.OrganismParts, Is.Empty);
            Assert.That(hit.YearlyDownloads, Is.Empty);
            Assert.That(hit.OtherOmicsLinks, Is.Empty);
            Assert.That(hit.Highlights, Is.Empty);
            Assert.That(hit.BotCount, Is.Zero);
            Assert.That(hit.AvgDownloadsPerFile, Is.Zero);
            Assert.That(hit.SubmissionDate, Is.EqualTo(default(DateTime)));
        });
    }

    // ---- paging -------------------------------------------------------------

    [Test]
    public async Task SearchProjectsAsync_PagesUntilTotalRecords_ConcatenatesAllHits()
    {
        const long total = 3;
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("page=0")) return JsonResponse(Array(HitJson("PXD000001"), HitJson("PXD000002")), totalRecords: total);
            if (uri.Contains("page=1")) return JsonResponse(Array(HitJson("PXD000003")), totalRecords: total);
            return JsonResponse("[]", totalRecords: total);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("liver", pageSize: 2);

        Assert.That(hits.Select(h => h.Accession), Is.EqualTo(new[] { "PXD000001", "PXD000002", "PXD000003" }));
        // Stopped at total_records on a visibly short final page: no tail probe, no wasted request.
        Assert.That(handler.RequestedUris.Count, Is.EqualTo(2));
    }

    /// <summary>
    /// Search inherits the shared pager, so it inherits the understated-total_records defence too:
    /// a fetch ending on a page that could be full is probed once rather than trusted. Pinned here
    /// because the fix was made for the manifest and nothing else would notice if search ever stopped
    /// consuming it.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_UnderstatedTotalRecords_StillReturnsEveryHit()
    {
        const long understated = 2;
        var handler = new StubHandler(request =>
        {
            string uri = request.RequestUri.ToString();
            if (uri.Contains("page=0")) return JsonResponse(Array(HitJson("PXD000001"), HitJson("PXD000002")), totalRecords: understated);
            if (uri.Contains("page=1")) return JsonResponse(Array(HitJson("PXD000003")), totalRecords: understated);
            return JsonResponse("[]", totalRecords: understated);
        });
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("liver", pageSize: 2);

        Assert.That(hits.Select(h => h.Accession), Is.EqualTo(new[] { "PXD000001", "PXD000002", "PXD000003" }));
    }

    // ---- absence is not an error --------------------------------------------

    /// <summary>
    /// Zero hits is normal absence, not a failure: PRIDE answers 200 with an empty array and
    /// total_records 0. This is why there is no TrySearchProjectsAsync — the Try/throwing pair on
    /// GetProjectAsync exists for a 404 that this endpoint never returns.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_NoMatches_ReturnsEmptyListNotNull()
    {
        var handler = new StubHandler(_ => JsonResponse("[]", totalRecords: 0));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("no-project-matches-this-keyword");

        Assert.That(hits, Is.Not.Null);
        Assert.That(hits, Is.Empty);
        Assert.That(handler.RequestedUris.Count, Is.EqualTo(1));
    }

    // ---- request construction -----------------------------------------------

    [Test]
    public async Task SearchProjectsAsync_SendsKeywordAndPageSizeAndPage()
    {
        var handler = new StubHandler(_ => JsonResponse("[]", totalRecords: 0));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        await client.SearchProjectsAsync("liver disease", pageSize: 25);

        string uri = handler.RequestedUris.Single();
        Assert.Multiple(() =>
        {
            Assert.That(uri, Does.StartWith("https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?"));
            // The space IS escaped to %20 on the wire; Uri.ToString() decodes it back for display, so
            // the decoded form is what this assertion can see. Asserting "%20" here would be asserting
            // a property of Uri.ToString() rather than of the request — the characters that actually
            // restructure a query are covered in the escaping test below, where %26/%3D do survive.
            Assert.That(uri, Does.Contain("keyword=liver disease"));
            Assert.That(uri, Does.Contain("pageSize=25"));
            Assert.That(uri, Does.Contain("page=0"));
        });
    }

    /// <summary>
    /// A keyword is free text, so the characters that matter are the ones that would restructure the
    /// query string. An unescaped "&amp;" splits the keyword into a second parameter and an unescaped
    /// "=" splits it into a name/value pair — either way PRIDE answers 200 with results for a
    /// DIFFERENT query, so the caller gets wrong data rather than an error. Uri.ToString() preserves
    /// %26 and %3D (unlike %20 and %2F, which it decodes), so both are assertable directly.
    /// </summary>
    [Test]
    public async Task SearchProjectsAsync_KeywordWithQuerySyntax_IsEscaped()
    {
        var handler = new StubHandler(_ => JsonResponse("[ ]", totalRecords: 0));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        await client.SearchProjectsAsync("a&pageSize=1#b");

        string uri = handler.RequestedUris.Single();
        Assert.Multiple(() =>
        {
            Assert.That(uri, Does.Contain("keyword=a%26pageSize%3D1%23b"));
            // The injected pageSize must not have displaced the real one.
            Assert.That(uri, Does.Contain("pageSize=100"));
            Assert.That(uri.Split("pageSize=").Length - 1, Is.EqualTo(1), "the keyword must not add a second pageSize parameter");
            // '#' is the dangerous one, and the reason this test carries it: unescaped it is a FRAGMENT
            // delimiter, so HttpClient sends nothing after it -- truncating the keyword AND deleting
            // pageSize and page from the request, which then silently fall back to server defaults.
            // Verified live 2026-08-21: the response grew from 12 KB to ~1 MB. Losing `page` would
            // defeat paging itself, which is the failure GetAllPagesAsync exists to prevent.
            Assert.That(uri, Does.Contain("page=0"), "the keyword must not be able to strip the paging parameters");
        });
    }

    /// <summary>
    /// An over-long keyword is refused here rather than sent, because PRIDE answers it with HTTP 500
    /// (observed at 2000 characters on 2026-08-21) rather than a 400 or a 414.
    /// </summary>
    /// <remarks>
    /// A 500 is the one status that cannot be told apart from an outage: ExternalServiceTestHelper
    /// reads it as "PRIDE is down" and SKIPS, and a caller reading the exception message sees the same
    /// thing. So without this guard a caller's own bug would masquerade as an EBI outage -- and a live
    /// test covering it would go green by skipping. Failing fast at the call site is the same posture
    /// the method already takes for a blank keyword and a non-positive page size.
    /// </remarks>
    [Test]
    public void SearchProjectsAsync_KeywordLongerThanTheCap_ThrowsWithoutSendingAnything()
    {
        var handler = new StubHandler(_ => JsonResponse("[ ]", totalRecords: 0));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        Assert.ThrowsAsync<ArgumentException>(
            () => client.SearchProjectsAsync(new string('a', PrideArchiveClient.MaxKeywordLength + 1)));
        Assert.That(handler.RequestedUris, Is.Empty, "the request must not reach PRIDE at all");
    }

    [Test]
    public async Task SearchProjectsAsync_KeywordExactlyAtTheCap_IsSent()
    {
        var handler = new StubHandler(_ => JsonResponse("[ ]", totalRecords: 0));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        await client.SearchProjectsAsync(new string('a', PrideArchiveClient.MaxKeywordLength));

        Assert.That(handler.RequestedUris.Count, Is.EqualTo(1), "the cap is inclusive");
    }

    // ---- argument and transport contracts -----------------------------------

    [TestCase(null)]
    [TestCase("")]
    [TestCase("   ")]
    public void SearchProjectsAsync_BlankKeyword_Throws(string keyword)
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));

        // Not a degenerate search: PRIDE reads an absent keyword as "browse the whole archive", which
        // is a different capability with a different cost, so it must be asked for deliberately.
        Assert.ThrowsAsync<ArgumentException>(() => client.SearchProjectsAsync(keyword));
    }

    [TestCase(0)]
    [TestCase(-1)]
    public void SearchProjectsAsync_NonPositivePageSize_Throws(int pageSize)
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));

        Assert.ThrowsAsync<ArgumentOutOfRangeException>(() => client.SearchProjectsAsync("liver", pageSize));
    }

    /// <summary>
    /// A non-success status is a transport failure and stays an HttpRequestException — the type
    /// ExternalServiceTestHelper converts to a SKIPPED live test. Contract breaks throw MzLibException
    /// instead, precisely so an outage and a broken promise are never confused.
    /// </summary>
    [Test]
    public void SearchProjectsAsync_NonSuccessStatus_ThrowsHttpRequestException()
    {
        var handler = new StubHandler(_ => JsonResponse("", HttpStatusCode.InternalServerError));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        Assert.ThrowsAsync<HttpRequestException>(() => client.SearchProjectsAsync("liver"));
    }

    [Test]
    public void SearchProjectsAsync_ServerRepeatsAPageWhileTotalSaysMoreRemains_ThrowsMzLibException()
    {
        // A server that ignores `page` while total_records promises more is a broken contract, not an
        // outage, so it must not surface as the exception type that means "PRIDE is down".
        var handler = new StubHandler(_ => JsonResponse(Array(HitJson("PXD000001")), totalRecords: 500));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        Assert.ThrowsAsync<MzLibException>(() => client.SearchProjectsAsync("liver", pageSize: 1));
    }

    [Test]
    public void SearchProjectsAsync_Cancelled_Throws()
    {
        var handler = new StubHandler(_ => JsonResponse(Array(HitJson("PXD000001")), totalRecords: 1));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.ThrowsAsync<OperationCanceledException>(() => client.SearchProjectsAsync("liver", 100, cts.Token));
    }
}

/// <summary>
/// Live canary against the real PRIDE Archive search endpoint. Carries [Category("ExternalService")]
/// so CI runs it in the dedicated, non-blocking external-service job rather than the required
/// unit-test run; [Category("Pride")] allows selecting it on its own. Routed through
/// <see cref="ExternalServiceTestHelper.RunAsync"/> so a PRIDE outage SKIPS the test while a genuine
/// contract break still FAILS.
/// </summary>
/// <remarks>
/// The drift worth catching here is the flattening. The offline fixture is a frozen copy of the wire
/// format and so cannot notice if PRIDE ever starts serving structured CvParam objects from this
/// endpoint — which would make every collection on the DTO deserialize to nothing. Only a live call
/// can, so that is what the assertions below are aimed at.
/// </remarks>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideProjectSearchLiveTests
{
    /// <summary>
    /// The keyword is deliberately a NARROW one. SearchProjectsAsync fetches every page, so a broad
    /// term ("liver" was 2 197 hits when this was written) would page through the whole result set on
    /// every CI run — hammering EBI to prove something a small result set proves just as well.
    /// "MetaMorpheus" was 93 hits, so this is a single request, and it is a term the lab's own output
    /// keeps alive.
    /// </summary>
    [Test]
    public Task SearchProjectsAsync_LiveKeyword_ReturnsPopulatedHits() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("MetaMorpheus");

            Assert.That(hits, Is.Not.Empty);
            Assert.Multiple(() =>
            {
                // Constraint form, not Assert.That(bool): a collapsed LINQ predicate reports only
                // "Expected: True, But was: False", which cannot be triaged from a CI log.
                Assert.That(hits, Is.All.Matches<PrideProjectSearchResult>(h => !string.IsNullOrEmpty(h.Accession)));
                Assert.That(hits, Is.All.Matches<PrideProjectSearchResult>(h => !string.IsNullOrEmpty(h.Title)));
                Assert.That(hits, Is.All.Matches<PrideProjectSearchResult>(h => h.SubmissionDate != default));

                // The flattening itself: if PRIDE ever returns objects here these come back empty.
                Assert.That(hits.Any(h => h.Organisms.Count > 0), Is.True, "no hit carried an organism");
                Assert.That(hits.Any(h => h.Instruments.Count > 0), Is.True, "no hit carried an instrument");
            });
        });

    [Test]
    public Task SearchProjectsAsync_LiveNoMatches_ReturnsEmpty() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            List<PrideProjectSearchResult> hits =
                await client.SearchProjectsAsync("zzzzqqqqxxxx-no-such-dataset-zzzzqqqqxxxx");

            // Absence, not an error — the fact that removes the need for a Try variant.
            Assert.That(hits, Is.Not.Null);
            Assert.That(hits, Is.Empty);
        });

    /// <summary>
    /// Pins the multi-page path against the real server — the only place the server-side pageSize cap
    /// and PRIDE's real total_records are in play. Same narrow keyword as above (93 hits when written)
    /// at a small page size, so this costs about five requests rather than paging a broad result set.
    /// </summary>
    [Test]
    public Task SearchProjectsAsync_LiveMultiPage_ReturnsMoreThanOnePageOfHits() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            List<PrideProjectSearchResult> hits = await client.SearchProjectsAsync("MetaMorpheus", pageSize: 20);

            Assert.That(hits.Count, Is.GreaterThan(20), "search stopped at the first page instead of paging");
            Assert.That(hits.Select(h => h.Accession).Distinct().Count(), Is.EqualTo(hits.Count), "the pager returned a duplicate hit");
        });
}
