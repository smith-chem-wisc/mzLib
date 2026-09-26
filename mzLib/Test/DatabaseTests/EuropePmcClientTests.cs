using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;
using NUnit.Framework;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    /// <summary>
    /// The Europe PMC client, offline against captured payloads (the PRIDE client's test convention), plus one live
    /// canary that skips on an outage.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class EuropePmcClientTests
    {
        private sealed class StubHandler : HttpMessageHandler
        {
            private readonly Func<HttpRequestMessage, HttpResponseMessage> _responder;
            public List<string> RequestedUris { get; } = new();

            public StubHandler(Func<HttpRequestMessage, HttpResponseMessage> responder) => _responder = responder;

            protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                RequestedUris.Add(request.RequestUri!.ToString());
                return Task.FromResult(_responder(request));
            }
        }

        private static HttpResponseMessage Text(string body, HttpStatusCode status = HttpStatusCode.OK, string type = "application/json") =>
            new(status) { Content = new StringContent(body, Encoding.UTF8, type) };

        private static HttpResponseMessage Bytes(byte[] body) =>
            new(HttpStatusCode.OK) { Content = new ByteArrayContent(body) };

        // Captured live on 2026-09-26 from search?query=EXT_ID:31836719 AND SRC:MED&resultType=lite, trimmed.
        private const string Hit = """
            {"version":"6.9","hitCount":1,"resultList":{"result":[{"id":"31836719","source":"MED","pmid":"31836719","pmcid":"PMC6910996",
            "fullTextIdList":{"fullTextId":["PMC6910996"]},"doi":"10.1038/s41597-019-0317-x",
            "title":"Engaging a community to enable disease-centric data sharing with the NF Data Portal.",
            "journalTitle":"Sci Data","pubYear":"2019","isOpenAccess":"Y","inEPMC":"Y","inPMC":"Y","hasSuppl":"N"}]}}
            """;

        // Captured live on 2026-09-26: a search with no match.
        private const string NoHit = """{"version":"6.9","hitCount":0,"resultList":{"result":[]}}""";

        // Captured live on 2026-09-26 from PMC3650345/supplementaryFiles: Europe PMC's "none" is a 200 with this body.
        private const string NoSupplement = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?><ns4:errorBean xmlns:ns4="http://webservice.cdb.ebi.ac.uk/"><errCode>0</errCode><errMsg>Article with id PMC3650345 is not open access one</errMsg></ns4:errorBean>""";

        private static string TempDir()
        {
            string d = Path.Combine(TestContext.CurrentContext.WorkDirectory, "EuropePmc_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(d);
            return d;
        }

        private static byte[] Zip()
        {
            using var ms = new MemoryStream();
            using (var zip = new ZipArchive(ms, ZipArchiveMode.Create, leaveOpen: true))
                using (var w = new StreamWriter(zip.CreateEntry("mmc1.csv").Open())) w.Write("Sample,Sex\nA,female\n");
            return ms.ToArray();
        }

        [Test]
        public async Task APubMedIdFindsTheArticleAndWhetherItsTextIsOpen()
        {
            var handler = new StubHandler(_ => Text(Hit));
            using var client = new EuropePmcClient(new HttpClient(handler));

            var (found, article) = await client.TryFindArticleAsync(pubMedId: 31836719, doi: "");

            Assert.That(found, Is.True);
            Assert.That((article.PubMedId, article.PmcId, article.Doi), Is.EqualTo(("31836719", "PMC6910996", "10.1038/s41597-019-0317-x")));
            Assert.That(article.IsOpenAccess, Is.True);
            Assert.That(article.HasFullText, Is.True);
            Assert.That(handler.RequestedUris.Single(), Does.Contain("EXT_ID%3A31836719"));
        }

        [Test]
        public async Task ADoiIsTriedWhenThereIsNoPubMedIdAndNoMatchIsAValue()
        {
            var handler = new StubHandler(_ => Text(NoHit));
            using var client = new EuropePmcClient(new HttpClient(handler));

            var (found, article) = await client.TryFindArticleAsync(pubMedId: 0, doi: "10.1097/HC9.0000000000000795");

            Assert.That(found, Is.False);
            Assert.That(article, Is.Null);
            Assert.That(Uri.UnescapeDataString(handler.RequestedUris.Single()), Does.Contain("DOI:\"10.1097/HC9.0000000000000795\""), "no PubMed id, so the DOI is searched");
        }

        [Test]
        public void NothingToSearchForIsAnArgumentErrorBeforeAnyRequest()
        {
            var handler = new StubHandler(_ => Text(Hit));
            using var client = new EuropePmcClient(new HttpClient(handler));

            Assert.ThrowsAsync<ArgumentException>(() => client.TryFindArticleAsync(0, " "));
            Assert.That(handler.RequestedUris, Is.Empty);
        }

        [TestCase(HttpStatusCode.ServiceUnavailable, typeof(HttpRequestException))]
        [TestCase(HttpStatusCode.TooManyRequests, typeof(HttpRequestException))]
        [TestCase(HttpStatusCode.Forbidden, typeof(MzLibException))]
        public void AnOutageIsATransportFailureAndAnyOtherRefusalIsNot(HttpStatusCode status, Type expected)
        {
            using var client = new EuropePmcClient(new HttpClient(new StubHandler(_ => Text("{}", status))));

            var e = Assert.CatchAsync(() => client.TryFindArticleAsync(31836719, ""));
            Assert.That(e, Is.TypeOf(expected));
            Assert.That(e!.Message, Does.Not.Contain("http"), "no URL in the message (mzLib #1350)");
        }

        [Test]
        public async Task TheFullTextIsSavedAndReusedAndNotRequestedWhenClosed()
        {
            string dir = TempDir();
            var handler = new StubHandler(_ => Text("<?xml version=\"1.0\"?><article><body><sec><title>Methods</title></sec></body></article>", type: "application/xml"));
            using var client = new EuropePmcClient(new HttpClient(handler));
            var article = new EuropePmcArticle { PubMedId = "31836719", PmcId = "PMC6910996", IsOpenAccess = true, HasFullText = true };

            var (found, path) = await client.TryDownloadFullTextXmlAsync(article, dir);
            var (again, samePath) = await client.TryDownloadFullTextXmlAsync(article, dir);

            Assert.That((found, again), Is.EqualTo((true, true)));
            Assert.That(Path.GetFileName(path), Is.EqualTo("PMC6910996.xml"));
            Assert.That(samePath, Is.EqualTo(path));
            Assert.That(handler.RequestedUris, Has.Count.EqualTo(1), "the second call reuses the file on disk");

            // Europe PMC answers a closed article's full text with HTTP 500, which would read as an outage: never ask.
            var closed = new EuropePmcArticle { PubMedId = "1", PmcId = "", HasFullText = false };
            Assert.That((await client.TryDownloadFullTextXmlAsync(closed, dir)).Found, Is.False);
            Assert.That(handler.RequestedUris, Has.Count.EqualTo(1));
        }

        [Test]
        public async Task SupplementsAreSavedOnlyWhenTheyAreAZip()
        {
            string dir = TempDir();
            var article = new EuropePmcArticle { PmcId = "PMC3650345", HasFullText = true };

            using (var none = new EuropePmcClient(new HttpClient(new StubHandler(_ => Text(NoSupplement, type: "application/xml")))))
            {
                var (found, _) = await none.TryDownloadSupplementaryFilesAsync(article, dir);
                Assert.That(found, Is.False);
                Assert.That(Directory.GetFiles(dir), Is.Empty, "an error page is never cached");
            }

            using var some = new EuropePmcClient(new HttpClient(new StubHandler(_ => Bytes(Zip()))));
            var (ok, path) = await some.TryDownloadSupplementaryFilesAsync(article, dir);
            Assert.That(ok, Is.True);
            Assert.That(Path.GetFileName(path), Is.EqualTo("PMC3650345_supplementary.zip"));
            using var zip = ZipFile.OpenRead(path);
            Assert.That(zip.Entries.Single().Name, Is.EqualTo("mmc1.csv"));
        }

        [Test]
        public async Task APrideReferenceIsResolvedThroughItsPubMedIdOrDoi()
        {
            var handler = new StubHandler(_ => Text(Hit));
            using var client = new EuropePmcClient(new HttpClient(handler));
            var reference = new PrideReference { PubmedId = 31836719, Doi = "10.1038/s41597-019-0317-x" };

            var (found, article) = await client.TryFindArticleAsync(reference);

            Assert.That(found, Is.True);
            Assert.That(article.PmcId, Is.EqualTo("PMC6910996"));
        }
    }

    [TestFixture]
    [Category("ExternalService")]
    [Category("EuropePmc")]
    [ExcludeFromCodeCoverage]
    public class EuropePmcClientLiveTests
    {
        [Test]
        public async Task ALivePubMedIdResolvesToAnOpenAccessArticle()
        {
            await ExternalServiceTestHelper.RunAsync("Europe PMC", async () =>
            {
                using var client = new EuropePmcClient();
                var (found, article) = await client.TryFindArticleAsync(31836719, "");
                Assert.That(found, Is.True);
                Assert.That(article.PmcId, Is.EqualTo("PMC6910996"));
                Assert.That(article.HasFullText, Is.True);
            });
        }
    }
}
