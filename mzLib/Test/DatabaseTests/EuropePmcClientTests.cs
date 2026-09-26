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

        /// <summary>Delivers a few bytes and then goes silent, honouring only cancellation (PrideArchiveDownloadTests' StallingStream).</summary>
        private sealed class StallingStream : Stream
        {
            private int _bytesBeforeStall;
            public StallingStream(int bytesBeforeStall) => _bytesBeforeStall = bytesBeforeStall;

            public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
            {
                if (_bytesBeforeStall > 0)
                {
                    int n = Math.Min(buffer.Length, _bytesBeforeStall);
                    buffer.Span.Slice(0, n).Fill((byte)'<');
                    _bytesBeforeStall -= n;
                    return n;
                }
                await Task.Delay(Timeout.Infinite, cancellationToken).ConfigureAwait(false);
                return 0;
            }

            public override int Read(byte[] buffer, int offset, int count) => ReadAsync(buffer.AsMemory(offset, count)).AsTask().GetAwaiter().GetResult();
            public override bool CanRead => true;
            public override bool CanSeek => false;
            public override bool CanWrite => false;
            public override long Length => throw new NotSupportedException();
            public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
            public override void Flush() { }
            public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
            public override void SetLength(long value) => throw new NotSupportedException();
            public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();
        }

        private static readonly EuropePmcArticle Open = new() { PmcId = "PMC6910996", HasFullText = true };

        [TestCase("<html>Service temporarily down</html>")]
        [TestCase("""{"version":"6.9","hitCount":0}""")]
        public void ASearchAnswerThatIsNotAResultListIsAContractFailure(string body)
        {
            using var client = new EuropePmcClient(new HttpClient(new StubHandler(_ => Text(body))));

            Assert.ThrowsAsync<MzLibException>(() => client.TryFindArticleAsync(31836719, ""));
        }

        [Test]
        public async Task AnArticleIsHeldInFullTextOnlyWhenEuropePmcSaysSo()
        {
            // Not in Europe PMC, with a PMCID Europe PMC does not list: no full text, so nothing is ever requested.
            const string notHeld = """{"resultList":{"result":[{"pmid":"1","pmcid":"PMC1","inEPMC":"N","fullTextIdList":{"fullTextId":["PMC9"]}}]}}""";
            // Listed without the inEPMC flag: held.
            const string listed = """{"resultList":{"result":[{"pmid":"2","pmcid":"PMC2","inEPMC":"N","fullTextIdList":{"fullTextId":["PMC2"]}}]}}""";
            // No identifiers at all: every string is empty, never null.
            const string bare = """{"resultList":{"result":[{"source":"MED"}]}}""";

            async Task<EuropePmcArticle> Find(string body)
            {
                using var client = new EuropePmcClient(new HttpClient(new StubHandler(_ => Text(body))));
                return (await client.TryFindArticleAsync(1, "")).Article;
            }

            Assert.That((await Find(notHeld)).HasFullText, Is.False);
            Assert.That((await Find(listed)).HasFullText, Is.True);
            var b = await Find(bare);
            Assert.That((b.PubMedId, b.PmcId, b.Doi, b.Title, b.HasFullText, b.IsOpenAccess), Is.EqualTo(("", "", "", "", false, false)));
        }

        [Test]
        public async Task AFullTextThatIsMissingIsAValueAndOneThatIsNotXmlIsAFailure()
        {
            string dir = TempDir();

            using (var missing = new EuropePmcClient(new HttpClient(new StubHandler(_ => Text("", HttpStatusCode.NotFound)))))
                Assert.That((await missing.TryDownloadFullTextXmlAsync(Open, dir)).Found, Is.False);

            foreach (var notXml in new[] { "{\"error\":true}", "   " })
            {
                using var client = new EuropePmcClient(new HttpClient(new StubHandler(_ => Text(notXml))));
                Assert.ThrowsAsync<MzLibException>(() => client.TryDownloadFullTextXmlAsync(Open, dir), notXml);
            }
            Assert.That(Directory.GetFiles(dir), Is.Empty, "nothing that failed is written");
        }

        [Test]
        public async Task AFullTextWithAByteOrderMarkIsXmlAndOverwriteFetchesItAgain()
        {
            string dir = TempDir();
            byte[] withBom = new byte[] { 0xEF, 0xBB, 0xBF }.Concat(Encoding.UTF8.GetBytes("\n <article/>")).ToArray();
            var handler = new StubHandler(_ => Bytes(withBom));
            using var client = new EuropePmcClient(new HttpClient(handler));

            var (found, path) = await client.TryDownloadFullTextXmlAsync(Open, dir);
            await client.TryDownloadFullTextXmlAsync(Open, dir, overwrite: true);

            Assert.That(found, Is.True);
            Assert.That(File.ReadAllBytes(path), Is.EqualTo(withBom));
            Assert.That(handler.RequestedUris, Has.Count.EqualTo(2), "overwrite ignores the file on disk");
        }

        [TestCase(null)]
        [TestCase(new byte[] { (byte)'P', (byte)'K', 3 })]
        [TestCase(new byte[] { (byte)'P', (byte)'X', 3, 4 })]
        [TestCase(new byte[] { (byte)'<', (byte)'K', 3, 4 })]
        public async Task OnlyAZipIsASupplement(byte[] body)
        {
            string dir = TempDir();
            using var client = new EuropePmcClient(new HttpClient(new StubHandler(_ => body == null ? Text("", HttpStatusCode.NotFound) : Bytes(body))));

            var (found, _) = await client.TryDownloadSupplementaryFilesAsync(Open, dir);

            Assert.That(found, Is.False);
            Assert.That(Directory.GetFiles(dir), Is.Empty);
        }

        [Test]
        public async Task SupplementsOnDiskAreReusedAndAClosedArticleIsNeverAsked()
        {
            string dir = TempDir();
            var handler = new StubHandler(_ => Bytes(Zip()));
            using var client = new EuropePmcClient(new HttpClient(handler));

            await client.TryDownloadSupplementaryFilesAsync(Open, dir);
            var (again, _) = await client.TryDownloadSupplementaryFilesAsync(Open, dir);
            var (closed, _) = await client.TryDownloadSupplementaryFilesAsync(new EuropePmcArticle { PmcId = "PMC1" }, dir);
            var (noId, _) = await client.TryDownloadSupplementaryFilesAsync(new EuropePmcArticle { HasFullText = true }, dir);

            Assert.That((again, closed, noId), Is.EqualTo((true, false, false)));
            Assert.That(handler.RequestedUris, Has.Count.EqualTo(1));
        }

        [TestCase("PMC1/../../x")]
        [TestCase("..")]
        public void APmcIdThatIsAPathIsRefusedBeforeAnyRequest(string pmcId)
        {
            var handler = new StubHandler(_ => Bytes(Zip()));
            using var client = new EuropePmcClient(new HttpClient(handler));
            var article = new EuropePmcArticle { PmcId = pmcId, HasFullText = true };

            Assert.ThrowsAsync<ArgumentException>(() => client.TryDownloadSupplementaryFilesAsync(article, TempDir()));
            Assert.That(handler.RequestedUris, Is.Empty);
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => new EuropePmcClient(null!));
            using var client = new EuropePmcClient(new HttpClient(new StubHandler(_ => Bytes(Zip()))));
            Assert.ThrowsAsync<ArgumentNullException>(() => client.TryFindArticleAsync((PrideReference)null!));
            Assert.ThrowsAsync<ArgumentNullException>(() => client.TryDownloadFullTextXmlAsync(null!, TempDir()));
            Assert.ThrowsAsync<ArgumentException>(() => client.TryDownloadFullTextXmlAsync(Open, " "));
        }

        [Test]
        public async Task ASuppliedBaseAddressIsKept()
        {
            var handler = new StubHandler(_ => Text(Hit));
            using var client = new EuropePmcClient(new HttpClient(handler) { BaseAddress = new Uri("http://mirror.example/rest/") });

            await client.TryFindArticleAsync(31836719, "");

            Assert.That(handler.RequestedUris.Single(), Does.StartWith("http://mirror.example/rest/search"));
        }

        [Test]
        public void ABodyThatStallsIsAnOutageAndLeavesNoFile()
        {
            string dir = TempDir();
            var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new StallingStream(4)) });
            using var client = new EuropePmcClient(new HttpClient(handler)) { BodyStallTimeout = TimeSpan.FromMilliseconds(200) };

            var e = Assert.ThrowsAsync<HttpRequestException>(() => client.TryDownloadFullTextXmlAsync(Open, dir));

            Assert.That(e!.Message, Does.Contain("delivered nothing"));
            Assert.That(Directory.GetFiles(dir), Is.Empty);
        }

        [Test]
        public void ACallerCancellingDuringTheBodyStaysACancellation()
        {
            var handler = new StubHandler(_ => new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new StallingStream(4)) });
            using var client = new EuropePmcClient(new HttpClient(handler)) { BodyStallTimeout = TimeSpan.FromSeconds(30) };
            using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(200));

            Assert.CatchAsync<OperationCanceledException>(() => client.TryDownloadFullTextXmlAsync(Open, TempDir(), cancellationToken: cts.Token));
        }

        [Test]
        public void TheDefaultClientDisposesOnceAndASuppliedHttpClientIsLeftOpen()
        {
            var own = new EuropePmcClient();
            Assert.That(own.BodyStallTimeout, Is.EqualTo(TimeSpan.FromMinutes(2)));
            own.Dispose();
            Assert.That(() => own.Dispose(), Throws.Nothing, "idempotent");

            var http = new HttpClient(new StubHandler(_ => Text(Hit)));
            new EuropePmcClient(http).Dispose();
            Assert.That(() => http.GetAsync("http://x/").GetAwaiter().GetResult(), Throws.Nothing, "the caller's HttpClient is not disposed");
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
