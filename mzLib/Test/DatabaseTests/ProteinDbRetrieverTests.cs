using MzLibUtil;
using NUnit.Framework;
using Proteomics;
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
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// Offline tests for <see cref="ProteinDbRetriever"/>'s failure contract. Every UniProt response these
/// need — an outage, a 404, a query that matched nothing, a deleted accession — is one a live test could
/// not summon on demand, so they are driven through a stub handler and run in the required CI job.
/// The live counterpart is <see cref="ProteinDbRetrieverLiveTests"/>.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class ProteinDbRetrieverTests
{
    private string _storageDirectory;

    [SetUp]
    public void SetUp()
    {
        _storageDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteinDbRetrieverTests_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_storageDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_storageDirectory))
            Directory.Delete(_storageDirectory, recursive: true);
    }

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

    /// <summary>Response content whose stream fails part-way through, standing in for a dropped connection.</summary>
    private sealed class ThrowingContent : HttpContent
    {
        protected override Task SerializeToStreamAsync(Stream stream, TransportContext context) =>
            throw new IOException("the connection was reset");

        // The download path reads the body synchronously (ReadAsStream), so the sync override is the one
        // that actually runs; without it the content throws NotSupportedException instead of failing the
        // way a dropped connection does.
        protected override void SerializeToStream(Stream stream, TransportContext context, CancellationToken cancellationToken) =>
            throw new IOException("the connection was reset");

        protected override bool TryComputeLength(out long length)
        {
            length = 0;
            return false;
        }
    }

    private static HttpResponseMessage Response(string body, HttpStatusCode status = HttpStatusCode.OK,
        long? totalResults = null, string mediaType = "text/plain")
    {
        var response = new HttpResponseMessage(status) { Content = new StringContent(body, Encoding.UTF8, mediaType) };
        if (totalResults.HasValue)
            response.Headers.TryAddWithoutValidation("X-Total-Results", totalResults.Value.ToString());
        return response;
    }

    private static HttpClient ClientReturning(string body, HttpStatusCode status = HttpStatusCode.OK, long? totalResults = null) =>
        new(new StubHandler(_ => Response(body, status, totalResults)));

    /// <summary>
    /// A client for the two-request proteome flow: the "size=0" count probe answers with an empty body and
    /// an X-Total-Results header, and the "stream" request answers with the payload. Anything that varies
    /// only one of the two is expressed by overriding that argument.
    /// </summary>
    private static StubHandler ProteomeHandler(long totalResults, string payload,
        HttpStatusCode countStatus = HttpStatusCode.OK, HttpStatusCode streamStatus = HttpStatusCode.OK,
        bool omitCountHeader = false) =>
        new(request => request.RequestUri.ToString().Contains("/search?")
            ? Response("", countStatus, omitCountHeader ? null : totalResults)
            : Response(payload, streamStatus));

    private static HttpClient ProteomeClient(long totalResults, string payload,
        HttpStatusCode countStatus = HttpStatusCode.OK, HttpStatusCode streamStatus = HttpStatusCode.OK,
        bool omitCountHeader = false) =>
        new(ProteomeHandler(totalResults, payload, countStatus, streamStatus, omitCountHeader));

    /// <summary>The nginx page the "/proteome/search" URL used to return, and which used to be saved as a proteome.</summary>
    private const string NotFoundHtml =
        "<html>\r\n<head><title>404 Not Found</title></head>\r\n<body>\r\n<center><h1>404 Not Found</h1></center>\r\n</body>\r\n</html>";

    /// <summary>The body UniProt returns for a deleted accession: the wrapper and its copyright, no entry.</summary>
    private const string EntrylessUniProtXml =
        """
        <?xml version="1.0" encoding="UTF-8"  standalone="no" ?>
        <uniprot xmlns="http://uniprot.org/uniprot">
        <copyright>
        Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
        </copyright>
        </uniprot>
        """;

    /// <summary>A one-entry UniProtKB XML in the real wire shape, small enough to inline.</summary>
    private const string OneEntryUniProtXml =
        """
        <?xml version="1.0" encoding="UTF-8"  standalone="no" ?>
        <uniprot xmlns="http://uniprot.org/uniprot">
        <entry dataset="Swiss-Prot">
          <accession>P02768</accession>
          <name>ALBU_HUMAN</name>
          <protein><recommendedName><fullName>Serum albumin</fullName></recommendedName></protein>
          <gene><name type="primary">ALB</name></gene>
          <organism><name type="scientific">Homo sapiens</name></organism>
          <sequence length="10" mass="1149">PEPTIDEMRK</sequence>
        </entry>
        </uniprot>
        """;

    private const string OneEntryUniProtFasta = ">sp|P02768|ALBU_HUMAN Serum albumin OS=Homo sapiens\nPEPTIDEMRK\n";

    // ---- RetrieveProteome: the call is wrong (no request is made) ------------

    [Test]
    [TestCase(null)]
    [TestCase("")]
    [TestCase("   ")]
    public void RetrieveProteome_BlankProteomeId_ThrowsArgumentException(string proteomeId)
    {
        var handler = new StubHandler(_ => Response("should not be requested"));
        var ex = Assert.Throws<ArgumentException>(() => ProteinDbRetriever.RetrieveProteome(proteomeId, _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));

        Assert.That(ex.ParamName, Is.EqualTo("proteomeID"));
        Assert.That(handler.RequestedUris, Is.Empty, "a bad call must not reach UniProt");
    }

    [Test]
    public void RetrieveProteome_BlankStorageDirectory_ThrowsArgumentException()
    {
        var ex = Assert.Throws<ArgumentException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", "  ",
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ClientReturning("unused")));

        Assert.That(ex.ParamName, Is.EqualTo("absolutePathToStorageDirectory"));
    }

    /// <summary>
    /// An unsupported format is a programming error, not an absence — the distinction that used to be lost
    /// because both this and a UniProt outage returned null.
    /// </summary>
    [Test]
    [TestCase(ProteinDbRetriever.ProteomeFormat.gff)]
    [TestCase(ProteinDbRetriever.ProteomeFormat.rdf)]
    [TestCase(ProteinDbRetriever.ProteomeFormat.html)]
    public void RetrieveProteome_UnsupportedFormat_ThrowsArgumentOutOfRange(ProteinDbRetriever.ProteomeFormat format)
    {
        var handler = new StubHandler(_ => Response("should not be requested"));
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            format, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));

        Assert.That(ex.ParamName, Is.EqualTo("format"));
        Assert.That(ex.Message, Does.Contain(format.ToString()));
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    /// <summary>A missing storage directory is a local configuration problem, distinct from every remote failure.</summary>
    [Test]
    public void RetrieveProteome_MissingStorageDirectory_ThrowsDirectoryNotFound()
    {
        string missing = Path.Combine(_storageDirectory, "no-such-subdirectory");
        var handler = new StubHandler(_ => Response("should not be requested"));

        Assert.Throws<DirectoryNotFoundException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", missing,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    // ---- RetrieveProteome: the service is having a bad day -------------------

    /// <summary>
    /// The statuses that mean "try again later". These are the ONLY ones a live test may skip on, so
    /// nothing else may be reported as HttpRequestException.
    /// </summary>
    [Test]
    [TestCase(HttpStatusCode.RequestTimeout)]      // 408
    [TestCase(HttpStatusCode.TooManyRequests)]     // 429
    [TestCase(HttpStatusCode.InternalServerError)] // 500
    [TestCase(HttpStatusCode.BadGateway)]          // 502
    [TestCase(HttpStatusCode.ServiceUnavailable)]  // 503
    [TestCase(HttpStatusCode.GatewayTimeout)]      // 504
    public void RetrieveProteome_ServiceUnavailable_ThrowsHttpRequestException(HttpStatusCode status)
    {
        var ex = Assert.Throws<HttpRequestException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ClientReturning("upstream is down", status)));

        Assert.That(ex.Message, Does.Contain(((int)status).ToString()));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty, "a failed retrieval must leave no file behind");
    }

    /// <summary>
    /// A client-side timeout arrives as TaskCanceledException; it means the same thing as a 504 and must be
    /// reported the same way, so a caller has one exception type to catch for "unavailable".
    /// </summary>
    [Test]
    public void RetrieveProteome_Timeout_ThrowsHttpRequestException()
    {
        var handler = new StubHandler(_ => throw new TaskCanceledException("timed out"));

        Assert.Throws<HttpRequestException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));
    }

    // ---- RetrieveProteome: UniProt answered, and there is nothing there ------

    /// <summary>
    /// The regression this whole change exists for: the XML branch pointed at a URL that no longer exists,
    /// so UniProt's 404 page was written to disk under a .xml name and its path returned as a success. It
    /// must now throw — and NOT as HttpRequestException, which a live test would skip and never notice.
    /// </summary>
    [Test]
    public void RetrieveProteome_NotFound_ThrowsMzLibExceptionAndWritesNothing()
    {
        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ClientReturning(NotFoundHtml, HttpStatusCode.NotFound)));

        Assert.That(ex.Message, Does.Contain("404"));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty, "the 404 page must not be saved as a proteome");
    }

    /// <summary>
    /// UniProt answers an unknown proteome ID with 200 and X-Total-Results: 0. Without this check the caller
    /// gets a path to an empty file that loads as zero proteins with nothing explaining why.
    /// </summary>
    [Test]
    public void RetrieveProteome_ZeroResults_ThrowsMzLibExceptionAndWritesNothing()
    {
        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveProteome("UP999999999", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ClientReturning("", HttpStatusCode.OK, totalResults: 0)));

        Assert.That(ex.Message, Does.Contain("UP999999999"));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    /// <summary>
    /// Without the count header there is no way to tell a real proteome from a mistyped one, and guessing
    /// "assume it is fine" is how an empty download came to be reported as a success. Say so instead.
    /// </summary>
    [Test]
    public void RetrieveProteome_NoTotalResultsHeader_ThrowsMzLibExceptionAndWritesNothing()
    {
        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ProteomeClient(4, OneEntryUniProtFasta, omitCountHeader: true)));

        Assert.That(ex.Message, Does.Contain("X-Total-Results"));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    /// <summary>
    /// The original bug, reachable through a different door: UniProt is known to answer a streaming request
    /// with 200 and an error message in the body. The count probe said there were proteins, so an empty
    /// download is a broken transfer — not something to write to disk and return the path of.
    /// </summary>
    [Test]
    public void RetrieveProteome_CountSaysNonEmptyButBodyIsEmpty_ThrowsMzLibExceptionAndWritesNothing()
    {
        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ProteomeClient(4, payload: "")));

        Assert.That(ex.Message, Does.Contain("empty body"));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty, "an empty download must not be left on disk");
    }

    /// <summary>
    /// A proteome ID becomes part of a file name. Without a guard, "../../evil" writes outside the storage
    /// directory, and an NTFS ':' silently creates an alternate data stream instead of a file — the same
    /// class of defect RetrieveEntry's accession guard exists to prevent.
    /// </summary>
    [Test]
    [TestCase("../../evil")]
    [TestCase(@"..\..\evil")]
    [TestCase("organism_id:9606")]
    public void RetrieveProteome_UnsafeProteomeId_ThrowsArgumentExceptionAndWritesNothing(string proteomeID)
    {
        var handler = ProteomeHandler(4, OneEntryUniProtFasta);

        var ex = Assert.Throws<ArgumentException>(() => ProteinDbRetriever.RetrieveProteome(proteomeID, _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));

        Assert.That(ex.ParamName, Is.EqualTo("proteomeID"));
        Assert.That(handler.RequestedUris, Is.Empty, "a bad call must not reach UniProt");
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    /// <summary>
    /// The interruption the ".partial" dance exists for. Without it the destination is opened directly, so a
    /// transfer that dies mid-body leaves a truncated file — or destroys a good one that was already there.
    /// </summary>
    [Test]
    public void RetrieveProteome_InterruptedDownload_LeavesTheExistingFileIntactAndNoPartial()
    {
        string destination = Path.Combine(_storageDirectory, "UP000008595_reviewed.fasta");
        File.WriteAllText(destination, "a previously good download");

        var handler = new StubHandler(request => request.RequestUri.ToString().Contains("/search?")
            ? Response("", HttpStatusCode.OK, totalResults: 4)
            : new HttpResponseMessage(HttpStatusCode.OK) { Content = new ThrowingContent() });

        Assert.Throws<HttpRequestException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));

        Assert.That(File.ReadAllText(destination), Is.EqualTo("a previously good download"),
            "an interrupted download must not damage the file already at the destination");
        Assert.That(Directory.GetFiles(_storageDirectory, "*.partial"), Is.Empty);
    }

    // ---- RetrieveProteome: success ------------------------------------------

    /// <summary>
    /// Pins BOTH URLs exactly, for every combination of the enums. Asserting only the endpoint path let the
    /// query semantics drift: the "_reviewed"/"_unreviewed" file name is built from a separate expression to
    /// the "reviewed:" query term, so a mutation making the query always-true would name the file
    /// "_unreviewed" and fill it with reviewed proteins without failing anything.
    /// </summary>
    [Test]
    // format,                                    reviewed,  compress, isoforms, expected count query,                        expected download query
    [TestCase(ProteinDbRetriever.ProteomeFormat.fasta, true,  false, true,
        "https://rest.uniprot.org/uniprotkb/search?query=UP000008595+AND+reviewed:true&format=list&size=0",
        "https://rest.uniprot.org/uniprotkb/stream?query=UP000008595+AND+reviewed:true&compressed=false&format=fasta&includeIsoforms:true")]
    [TestCase(ProteinDbRetriever.ProteomeFormat.fasta, false, true,  false,
        "https://rest.uniprot.org/uniprotkb/search?query=UP000008595+AND+reviewed:false&format=list&size=0",
        "https://rest.uniprot.org/uniprotkb/stream?query=UP000008595+AND+reviewed:false&compressed=true&format=fasta&includeIsoforms:false")]
    [TestCase(ProteinDbRetriever.ProteomeFormat.xml,   true,  true,  false,
        "https://rest.uniprot.org/uniprotkb/search?query=UP000008595+AND+reviewed:true&format=list&size=0",
        "https://rest.uniprot.org/uniprotkb/stream?query=UP000008595+AND+reviewed:true&compressed=true&format=xml")]
    [TestCase(ProteinDbRetriever.ProteomeFormat.xml,   false, false, true,
        "https://rest.uniprot.org/uniprotkb/search?query=UP000008595+AND+reviewed:false&format=list&size=0",
        "https://rest.uniprot.org/uniprotkb/stream?query=UP000008595+AND+reviewed:false&compressed=false&format=xml")]
    public void RetrieveProteome_CountsThenStreams_WithTheExactQuery(ProteinDbRetriever.ProteomeFormat format,
        bool reviewed, bool compress, bool isoforms, string expectedCountUrl, string expectedDownloadUrl)
    {
        var handler = ProteomeHandler(4, OneEntryUniProtFasta);

        ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory, format,
            reviewed ? ProteinDbRetriever.Reviewed.yes : ProteinDbRetriever.Reviewed.no,
            compress ? ProteinDbRetriever.Compress.yes : ProteinDbRetriever.Compress.no,
            isoforms ? ProteinDbRetriever.IncludeIsoforms.yes : ProteinDbRetriever.IncludeIsoforms.no,
            new HttpClient(handler));

        Assert.That(handler.RequestedUris, Has.Count.EqualTo(2), "one count probe, then one download");
        Assert.That(handler.RequestedUris[0], Is.EqualTo(expectedCountUrl));
        Assert.That(handler.RequestedUris[1], Is.EqualTo(expectedDownloadUrl));
    }

    /// <summary>
    /// The download must use "stream", never "search": "search" is paged and answers with only its first 25
    /// entries plus a cursor, so a human proteome came back as 25 of 20416 with nothing reporting a problem.
    /// </summary>
    [Test]
    public void RetrieveProteome_DownloadsFromStreamNotThePagedSearchEndpoint()
    {
        var handler = ProteomeHandler(20416, OneEntryUniProtFasta);

        ProteinDbRetriever.RetrieveProteome("UP000005640", _storageDirectory, ProteinDbRetriever.ProteomeFormat.fasta,
            ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no,
            new HttpClient(handler));

        Assert.That(handler.RequestedUris.Last(), Does.Contain("/uniprotkb/stream?"));
        Assert.That(handler.RequestedUris.Last(), Does.Not.Contain("/search?"));
    }

    /// <summary>A mistyped proteome ID must cost nothing: the count probe answers, and no download follows.</summary>
    [Test]
    public void RetrieveProteome_UnknownProteome_DoesNotAttemptTheDownload()
    {
        var handler = ProteomeHandler(0, OneEntryUniProtFasta);

        Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveProteome("UP999999999", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, new HttpClient(handler)));

        Assert.That(handler.RequestedUris, Has.Count.EqualTo(1));
        Assert.That(handler.RequestedUris.Single(), Does.Contain("size=0"));
    }

    /// <summary>The file names callers depend on, unchanged by this rework.</summary>
    [Test]
    [TestCase(ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.yes, ProteinDbRetriever.IncludeIsoforms.yes, "UP000008595_reviewed_isoform.fasta.gz")]
    [TestCase(ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.no, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no, "UP000008595_unreviewed.fasta")]
    [TestCase(ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.yes, ProteinDbRetriever.IncludeIsoforms.no, "UP000008595_reviewed.xml.gz")]
    [TestCase(ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.no, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no, "UP000008595_unreviewed.xml")]
    // xml carries no isoforms, so asking for them must not change the name
    [TestCase(ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.yes, "UP000008595_reviewed.xml")]
    public void RetrieveProteome_NamesTheFileAsBefore(ProteinDbRetriever.ProteomeFormat format,
        ProteinDbRetriever.Reviewed reviewed, ProteinDbRetriever.Compress compress,
        ProteinDbRetriever.IncludeIsoforms include, string expectedFileName)
    {
        string path = ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory, format, reviewed, compress,
            include, ClientReturning("payload", totalResults: 4));

        Assert.That(path, Is.EqualTo(Path.Combine(_storageDirectory, expectedFileName)));
        Assert.That(File.Exists(path));
    }

    /// <summary>
    /// Re-retrieving used to fail with an IOException, because the download opened the destination with
    /// FileMode.CreateNew. It now overwrites, and leaves no ".partial" scratch file behind either way.
    /// </summary>
    [Test]
    public void RetrieveProteome_OverwritesAnExistingFile_AndLeavesNoPartial()
    {
        string first = ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ProteomeClient(4, "first payload"));

        string second = ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.yes, ProteinDbRetriever.Compress.no,
            ProteinDbRetriever.IncludeIsoforms.no, ProteomeClient(4, "second payload"));

        Assert.That(second, Is.EqualTo(first));
        Assert.That(File.ReadAllText(second), Is.EqualTo("second payload"));
        Assert.That(Directory.GetFiles(_storageDirectory, "*.partial"), Is.Empty);
    }

    // ---- RetrieveEntry: the call is wrong ------------------------------------

    [Test]
    [TestCase(null)]
    [TestCase("")]
    [TestCase("   ")]
    public void RetrieveEntry_BlankAccession_ThrowsArgumentException(string accession)
    {
        var ex = Assert.Throws<ArgumentException>(() =>
            ProteinDbRetriever.RetrieveEntry(accession, _storageDirectory, ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("unused")));

        Assert.That(ex.ParamName, Is.EqualTo("accession"));
    }

    /// <summary>
    /// The accession becomes part of a file name, so one carrying path separators or ".." must be refused
    /// before it can write outside the storage directory.
    /// </summary>
    [Test]
    [TestCase("../../evil")]
    [TestCase(@"..\..\evil")]
    [TestCase("sub/P02768")]
    [TestCase("C:P02768")]
    [TestCase("P02768 ")]
    public void RetrieveEntry_UnsafeAccession_ThrowsArgumentExceptionAndWritesNothing(string accession)
    {
        var handler = new StubHandler(_ => Response(OneEntryUniProtXml));

        var ex = Assert.Throws<ArgumentException>(() =>
            ProteinDbRetriever.RetrieveEntry(accession, _storageDirectory, ProteinDbRetriever.ProteomeFormat.xml, new HttpClient(handler)));

        Assert.That(ex.ParamName, Is.EqualTo("accession"));
        Assert.That(handler.RequestedUris, Is.Empty);
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    [TestCase(ProteinDbRetriever.ProteomeFormat.gff)]
    [TestCase(ProteinDbRetriever.ProteomeFormat.txt)]
    public void RetrieveEntry_UnsupportedFormat_ThrowsArgumentOutOfRange(ProteinDbRetriever.ProteomeFormat format)
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() =>
            ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory, format, ClientReturning("unused")));

        Assert.That(ex.ParamName, Is.EqualTo("format"));
    }

    [Test]
    public void RetrieveEntry_MissingStorageDirectory_ThrowsDirectoryNotFound()
    {
        Assert.Throws<DirectoryNotFoundException>(() => ProteinDbRetriever.RetrieveEntry("P02768",
            Path.Combine(_storageDirectory, "no-such-subdirectory"), ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("unused")));
    }

    // ---- RetrieveEntry: no such entry, versus a bad day ----------------------

    /// <summary>
    /// UniProt answers a MALFORMED accession with 400 (not 404) and a well-formed unknown one with 404.
    /// Both mean "no such entry" and must be permanent — never HttpRequestException, which would be read as
    /// an outage and skip the live test that would otherwise catch a withdrawn accession.
    /// </summary>
    [Test]
    [TestCase(HttpStatusCode.BadRequest)] // 400: "NOTANACC"
    [TestCase(HttpStatusCode.NotFound)]   // 404: "P0ZZZ9"
    public void RetrieveEntry_NoSuchEntry_ThrowsMzLibException(HttpStatusCode status)
    {
        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveEntry("P0ZZZ9", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("<errorInfo/>", status)));

        Assert.That(ex.Message, Does.Contain("P0ZZZ9"));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    [TestCase(HttpStatusCode.RequestTimeout)]
    [TestCase(HttpStatusCode.TooManyRequests)]
    [TestCase(HttpStatusCode.ServiceUnavailable)]
    public void RetrieveEntry_ServiceUnavailable_ThrowsHttpRequestException(HttpStatusCode status)
    {
        Assert.Throws<HttpRequestException>(() => ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("down", status)));
    }

    /// <summary>
    /// A deleted accession is answered with 200 and a payload holding no entry at all. Saved as-is it loads
    /// as zero proteins, so it must be refused rather than reported as a successful retrieval.
    /// </summary>
    [Test]
    public void RetrieveEntry_TwoHundredWithNoEntry_ThrowsMzLibExceptionAndWritesNothing()
    {
        var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveEntry("Q9ZZZ9", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning(EntrylessUniProtXml)));

        Assert.That(ex.Message, Does.Contain("Q9ZZZ9"));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    public void RetrieveEntry_TwoHundredWithEmptyFasta_ThrowsMzLibException()
    {
        Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveEntry("Q9ZZZ9", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ClientReturning("")));

        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    /// <summary>
    /// 400 and 404 are the two statuses that mean "no such entry". Any OTHER non-success status is neither
    /// an absence nor an outage — ThrowIfServiceUnavailable has already claimed every transient status — so
    /// it must be MzLibException. Reporting it as HttpRequestException would make a live canary SKIP on a
    /// permanent contract break, which is the one outcome this whole design exists to prevent.
    /// </summary>
    [Test]
    [TestCase(HttpStatusCode.Unauthorized)]
    [TestCase(HttpStatusCode.Forbidden)]
    [TestCase(HttpStatusCode.UnsupportedMediaType)]
    public void RetrieveEntry_OtherNonSuccessStatus_ThrowsMzLibExceptionNotAnOutage(HttpStatusCode status)
    {
        Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("denied", status)));
    }

    // ---- TryRetrieveEntry: absence as a value, everything else still thrown ----

    /// <summary>
    /// The sibling of <see cref="ProteinDbRetriever.RetrieveEntry(string, string, ProteinDbRetriever.ProteomeFormat)"/>
    /// for an accession that came from user input, matching PrideArchiveClient's Get/TryGet pair. All three
    /// ways UniProt says "no such entry" become Found=false rather than an exception.
    /// </summary>
    [Test]
    [TestCase(HttpStatusCode.BadRequest, "<errorInfo/>")]                // 400: malformed accession
    [TestCase(HttpStatusCode.NotFound, "<errorInfo/>")]                  // 404: well-formed, unknown
    public void TryRetrieveEntry_NoSuchEntry_ReportsNotFoundWithoutThrowing(HttpStatusCode status, string body)
    {
        (bool found, string filePath) = ProteinDbRetriever.TryRetrieveEntry("P0ZZZ9", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning(body, status));

        Assert.That(found, Is.False);
        Assert.That(filePath, Is.Null);
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    public void TryRetrieveEntry_DeletedAccession_ReportsNotFoundWithoutThrowing()
    {
        (bool found, string filePath) = ProteinDbRetriever.TryRetrieveEntry("Q9ZZZ9", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning(EntrylessUniProtXml));

        Assert.That(found, Is.False);
        Assert.That(filePath, Is.Null);
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    public void TryRetrieveEntry_Success_ReportsFoundAndWritesTheFile()
    {
        (bool found, string filePath) = ProteinDbRetriever.TryRetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning(OneEntryUniProtXml));

        Assert.That(found, Is.True);
        Assert.That(filePath, Is.EqualTo(Path.Combine(_storageDirectory, "P02768.xml")));
        Assert.That(File.Exists(filePath));
    }

    /// <summary>
    /// An outage must NOT be collapsed into "no such entry" — that is the whole reason the Try form reports
    /// only one thing as a value.
    /// </summary>
    [Test]
    [TestCase(HttpStatusCode.ServiceUnavailable)]
    [TestCase(HttpStatusCode.TooManyRequests)]
    [TestCase(HttpStatusCode.RequestTimeout)]
    public void TryRetrieveEntry_ServiceUnavailable_StillThrows(HttpStatusCode status)
    {
        Assert.Throws<HttpRequestException>(() => ProteinDbRetriever.TryRetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("down", status)));
    }

    [Test]
    public void TryRetrieveEntry_BadArguments_StillThrow()
    {
        Assert.Throws<ArgumentException>(() => ProteinDbRetriever.TryRetrieveEntry("  ", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("unused")));
        Assert.Throws<ArgumentException>(() => ProteinDbRetriever.TryRetrieveEntry("../evil", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("unused")));
        Assert.Throws<DirectoryNotFoundException>(() => ProteinDbRetriever.TryRetrieveEntry("P02768",
            Path.Combine(_storageDirectory, "nope"), ProteinDbRetriever.ProteomeFormat.xml, ClientReturning("unused")));
    }

    // ---- RetrieveEntry: success ---------------------------------------------

    [Test]
    public void RetrieveEntry_RequestsTheEntryEndpointAndNamesTheFileForTheAccession()
    {
        var handler = new StubHandler(_ => Response(OneEntryUniProtXml, mediaType: "application/xml"));

        string path = ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, new HttpClient(handler));

        Assert.That(handler.RequestedUris.Single(), Is.EqualTo("https://rest.uniprot.org/uniprotkb/P02768.xml"));
        Assert.That(path, Is.EqualTo(Path.Combine(_storageDirectory, "P02768.xml")));
    }

    /// <summary>An isoform accession is a legitimate accession, and a legitimate file name.</summary>
    [Test]
    public void RetrieveEntry_IsoformAccession_IsAccepted()
    {
        var handler = new StubHandler(_ => Response(OneEntryUniProtXml, mediaType: "application/xml"));

        string path = ProteinDbRetriever.RetrieveEntry("P02768-2", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, new HttpClient(handler));

        Assert.That(handler.RequestedUris.Single(), Is.EqualTo("https://rest.uniprot.org/uniprotkb/P02768-2.xml"));
        Assert.That(path, Is.EqualTo(Path.Combine(_storageDirectory, "P02768-2.xml")));
    }

    /// <summary>
    /// The point of the method: what it writes is loadable by <see cref="ProteinDbLoader.LoadProteinXML"/>,
    /// which is what made reimplementing it downstream tempting in the first place.
    /// </summary>
    [Test]
    public void RetrieveEntry_WritesXmlThatProteinDbLoaderCanRead()
    {
        string path = ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.xml, ClientReturning(OneEntryUniProtXml));

        List<Protein> proteins = ProteinDbLoader.LoadProteinXML(path, generateTargets: true, DecoyType.None,
            null, isContaminant: false, null, out _);

        Assert.That(proteins.Count, Is.EqualTo(1));
        Assert.That(proteins[0].Accession, Is.EqualTo("P02768"));
        Assert.That(proteins[0].BaseSequence, Is.EqualTo("PEPTIDEMRK"));
    }

    [Test]
    public void RetrieveEntry_FastaFormat_WritesFastaThatProteinDbLoaderCanRead()
    {
        string path = ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory,
            ProteinDbRetriever.ProteomeFormat.fasta, ClientReturning(OneEntryUniProtFasta));

        Assert.That(path, Is.EqualTo(Path.Combine(_storageDirectory, "P02768.fasta")));

        List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(path, generateTargets: true, DecoyType.None,
            isContaminant: false, out _);

        Assert.That(proteins.Count, Is.EqualTo(1));
        Assert.That(proteins[0].BaseSequence, Is.EqualTo("PEPTIDEMRK"));
    }

    // ---- DownloadAvailableUniProtProteomes ----------------------------------

    [Test]
    public void DownloadAvailableUniProtProteomes_BlankFolder_ThrowsArgumentException()
    {
        var ex = Assert.Throws<ArgumentException>(() => ProteinDbRetriever.DownloadAvailableUniProtProteomes("  ", ClientReturning("unused")));
        Assert.That(ex.ParamName, Is.EqualTo("destinationFolder"));
    }

    /// <summary>
    /// "bubba" used to be indistinguishable from a UniProt outage: both returned null. It is a local
    /// configuration problem and now says so.
    /// </summary>
    [Test]
    public void DownloadAvailableUniProtProteomes_MissingFolder_ThrowsDirectoryNotFound()
    {
        Assert.Throws<DirectoryNotFoundException>(() =>
            ProteinDbRetriever.DownloadAvailableUniProtProteomes("bubba", ClientReturning("unused")));
    }

    [Test]
    [TestCase(HttpStatusCode.ServiceUnavailable)]
    [TestCase(HttpStatusCode.TooManyRequests)]
    public void DownloadAvailableUniProtProteomes_ServiceUnavailable_ThrowsHttpRequestException(HttpStatusCode status)
    {
        Assert.Throws<HttpRequestException>(() =>
            ProteinDbRetriever.DownloadAvailableUniProtProteomes(_storageDirectory, ClientReturning("down", status)));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    public void DownloadAvailableUniProtProteomes_NotFound_ThrowsMzLibException()
    {
        Assert.Throws<MzLibException>(() => ProteinDbRetriever.DownloadAvailableUniProtProteomes(
            _storageDirectory, ClientReturning(NotFoundHtml, HttpStatusCode.NotFound)));
        Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
    }

    [Test]
    public void DownloadAvailableUniProtProteomes_EmptyCatalogue_ThrowsMzLibException()
    {
        Assert.Throws<MzLibException>(() => ProteinDbRetriever.DownloadAvailableUniProtProteomes(
            _storageDirectory, ClientReturning("", HttpStatusCode.OK, totalResults: 0)));
    }

    [Test]
    public void DownloadAvailableUniProtProteomes_Success_WritesTheCatalogue()
    {
        string path = ProteinDbRetriever.DownloadAvailableUniProtProteomes(_storageDirectory,
            ClientReturning("gzip bytes would go here", totalResults: 2));

        Assert.That(path, Is.EqualTo(Path.Combine(_storageDirectory, "availableUniProtProteomes.txt.gz")));
        Assert.That(File.ReadAllText(path), Is.EqualTo("gzip bytes would go here"));
    }

    // ---- UniprotProteomesList -----------------------------------------------

    [Test]
    [TestCase(null)]
    [TestCase("  ")]
    public void UniprotProteomesList_BlankPath_ThrowsArgumentException(string path)
    {
        Assert.Throws<ArgumentException>(() => ProteinDbRetriever.UniprotProteomesList(path));
    }

    [Test]
    public void UniprotProteomesList_MissingFile_ThrowsFileNotFound()
    {
        var ex = Assert.Throws<FileNotFoundException>(() => ProteinDbRetriever.UniprotProteomesList("badFilePath"));
        Assert.That(ex.FileName, Is.EqualTo("badFilePath"));
    }

    /// <summary>An extension the reader cannot open is a bad argument, not an empty catalogue.</summary>
    [Test]
    public void UniprotProteomesList_UnsupportedExtension_ThrowsArgumentException()
    {
        string path = Path.Combine(_storageDirectory, "proteomes.fasta");
        File.WriteAllText(path, "UP000005640\tHomo sapiens (Human)");

        var ex = Assert.Throws<ArgumentException>(() => ProteinDbRetriever.UniprotProteomesList(path));
        Assert.That(ex.Message, Does.Contain(".fasta"));
    }

    [Test]
    public void UniprotProteomesList_ReadsGzippedTsv()
    {
        string path = Path.Combine(_storageDirectory, "availableUniProtProteomes.txt.gz");
        using (var fileStream = File.Create(path))
        using (var gzip = new GZipStream(fileStream, CompressionMode.Compress))
        using (var writer = new StreamWriter(gzip))
        {
            writer.WriteLine("Proteome Id\tOrganism");
            writer.WriteLine("UP000005640\tHomo sapiens (Human)");
        }

        Dictionary<string, string> proteomes = ProteinDbRetriever.UniprotProteomesList(path);

        Assert.That(proteomes["UP000005640"], Is.EqualTo("Homo sapiens (Human)"));
    }

    // ---- the classifiers themselves -----------------------------------------

    /// <summary>
    /// The boundary that everything else rests on: 4xx statuses that name a caller mistake must NOT be
    /// classified as an outage, or a live test would skip past a real regression.
    /// </summary>
    [Test]
    [TestCase(HttpStatusCode.OK, false)]
    [TestCase(HttpStatusCode.BadRequest, false)]
    [TestCase(HttpStatusCode.Unauthorized, false)]
    [TestCase(HttpStatusCode.Forbidden, false)]
    [TestCase(HttpStatusCode.NotFound, false)]
    [TestCase(HttpStatusCode.RequestTimeout, true)]
    [TestCase(HttpStatusCode.TooManyRequests, true)]
    [TestCase(HttpStatusCode.InternalServerError, true)]
    [TestCase(HttpStatusCode.ServiceUnavailable, true)]
    public void ThrowIfServiceUnavailable_ThrowsOnlyForServiceSideFailures(HttpStatusCode status, bool expectedThrow)
    {
        using var response = new HttpResponseMessage(status);

        if (expectedThrow)
            Assert.Throws<HttpRequestException>(() => ProteinDbRetriever.ThrowIfServiceUnavailable(response, "https://example.org"));
        else
            Assert.DoesNotThrow(() => ProteinDbRetriever.ThrowIfServiceUnavailable(response, "https://example.org"));
    }

    [Test]
    public void TryGetTotalResults_ReadsTheHeader()
    {
        using var response = new HttpResponseMessage(HttpStatusCode.OK);
        response.Headers.TryAddWithoutValidation("X-Total-Results", "4403");

        Assert.That(ProteinDbRetriever.TryGetTotalResults(response, out long total), Is.True);
        Assert.That(total, Is.EqualTo(4403));
    }

    [Test]
    [TestCase(null)]
    [TestCase("not a number")]
    public void TryGetTotalResults_AbsentOrUnparseableHeader_ReportsUnknown(string headerValue)
    {
        using var response = new HttpResponseMessage(HttpStatusCode.OK);
        if (headerValue != null)
            response.Headers.TryAddWithoutValidation("X-Total-Results", headerValue);

        Assert.That(ProteinDbRetriever.TryGetTotalResults(response, out _), Is.False);
    }
}

/// <summary>
/// Live canary against the real UniProt REST API. Carries [Category("ExternalService")] so CI runs it in
/// the dedicated, non-blocking external-service job (see .github/workflows/dotnet.yml) rather than the
/// required unit-test run; [Category("UniProt")] allows selecting it on its own. Calls are routed through
/// <see cref="ExternalServiceTestHelper.RunAsync"/>, so a UniProt outage (timeout / 5xx / unreachable)
/// SKIPS the test while a genuine contract break — a URL that has started 404ing, a response that no
/// longer parses, an accession that has been withdrawn — still FAILS. That separation is only possible
/// because <see cref="ProteinDbRetriever"/> now reports the two differently.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("UniProt")]
[ExcludeFromCodeCoverage]
public class ProteinDbRetrieverLiveTests
{
    private string _storageDirectory;

    [SetUp]
    public void SetUp()
    {
        _storageDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "ProteinDbRetrieverLive_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_storageDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_storageDirectory))
            Directory.Delete(_storageDirectory, recursive: true);
    }

    /// <summary>
    /// P02768 is human serum albumin: 609 residues, reviewed, and not going anywhere. The loader expands the
    /// entry's annotated sequence variants, so the assertions are on the canonical protein rather than on a
    /// count that grows whenever UniProt curates another variant.
    /// </summary>
    [Test]
    public Task RetrieveEntry_LiveP02768_LoadsAsOneAnnotatedProtein() =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            string path = ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory);

            Assert.That(path, Is.EqualTo(Path.Combine(_storageDirectory, "P02768.xml")));
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(path, generateTargets: true, DecoyType.None,
                null, isContaminant: false, null, out _);

            Assert.That(proteins, Is.Not.Empty);
            Assert.That(proteins.All(p => p.Accession.StartsWith("P02768")), "the file must hold one entry only");

            Protein canonical = proteins.Single(p => p.Accession == "P02768");
            Assert.That(canonical.BaseSequence.Length, Is.EqualTo(609));
            Assert.That(canonical.GeneNames.Any(g => g.Item2 == "ALB"));
            // The point of asking for xml rather than fasta: the annotations come with it.
            Assert.That(canonical.DatabaseReferences, Is.Not.Empty, "the entry should carry its cross-references");
            Assert.That(canonical.SequenceVariations, Is.Not.Empty, "the entry should carry its sequence variants");
            return Task.CompletedTask;
        });

    [Test]
    public Task RetrieveEntry_LiveFastaP02768_LoadsAsOneProtein() =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            string path = ProteinDbRetriever.RetrieveEntry("P02768", _storageDirectory, ProteinDbRetriever.ProteomeFormat.fasta);

            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(path, generateTargets: true, DecoyType.None,
                isContaminant: false, out _);

            Assert.That(proteins.Count, Is.EqualTo(1));
            Assert.That(proteins[0].BaseSequence.Length, Is.EqualTo(609));
            return Task.CompletedTask;
        });

    /// <summary>
    /// "NOTANACC" is malformed, so UniProt answers 400; "P0ZZZ9" is well-formed but unknown, so it answers
    /// 404. Both must arrive as a permanent MzLibException — if either came back as HttpRequestException the
    /// helper would skip this test and the distinction would go untested precisely when it mattered.
    /// </summary>
    [Test]
    [TestCase("NOTANACC")]
    [TestCase("P0ZZZ9")]
    public Task RetrieveEntry_LiveUnknownAccession_ThrowsMzLibException(string accession) =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveEntry(accession, _storageDirectory));
            Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
            return Task.CompletedTask;
        });

    /// <summary>
    /// The same three real absences through the Try form, which must report them as a value. Q9ZZZ9 is a
    /// deleted accession — UniProt answers 200 with a payload holding no entry.
    /// </summary>
    [Test]
    [TestCase("NOTANACC")] // 400
    [TestCase("P0ZZZ9")]   // 404
    [TestCase("Q9ZZZ9")]   // 200, no entry
    public Task TryRetrieveEntry_LiveUnknownAccession_ReportsNotFound(string accession) =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            (bool found, string filePath) = ProteinDbRetriever.TryRetrieveEntry(accession, _storageDirectory);

            Assert.That(found, Is.False);
            Assert.That(filePath, Is.Null);
            Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
            return Task.CompletedTask;
        });

    [Test]
    public Task TryRetrieveEntry_LiveP02768_ReportsFound() =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            (bool found, string filePath) = ProteinDbRetriever.TryRetrieveEntry("P02768", _storageDirectory);

            Assert.That(found, Is.True);
            Assert.That(filePath, Is.EqualTo(Path.Combine(_storageDirectory, "P02768.xml")));
            return Task.CompletedTask;
        });

    /// <summary>
    /// The XML branch of RetrieveProteome used to save a 404 page under this name. UP000008595 is Uukuniemi
    /// virus (strain S23), which has exactly 4 reviewed proteins.
    /// </summary>
    [Test]
    public Task RetrieveProteome_LiveXml_LoadsTheProteinsInsteadOfAnErrorPage() =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            string path = ProteinDbRetriever.RetrieveProteome("UP000008595", _storageDirectory,
                ProteinDbRetriever.ProteomeFormat.xml, ProteinDbRetriever.Reviewed.yes,
                ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no);

            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(path, generateTargets: true, DecoyType.None,
                null, isContaminant: false, null, out _);

            // Not an exact count: LoadProteinXML expands annotated sequence variants into extra Protein
            // objects, so curating one variant into any of the four entries would break an == 4 assertion
            // for a non-code reason. The four accessions are the stable fact.
            Assert.That(proteins.Select(p => p.Accession),
                Is.SupersetOf(new[] { "P09613", "P33453", "P22025", "P22026" }));
            return Task.CompletedTask;
        });

    /// <summary>
    /// Uukuniemi virus is fully reviewed, so asking for its unreviewed entries legitimately matches nothing.
    /// That used to be a zero-byte file whose path was returned as a success.
    /// </summary>
    [Test]
    public Task RetrieveProteome_LiveQueryMatchingNothing_ThrowsMzLibException() =>
        ExternalServiceTestHelper.RunAsync("UniProt", () =>
        {
            var ex = Assert.Throws<MzLibException>(() => ProteinDbRetriever.RetrieveProteome("UP000008595",
                _storageDirectory, ProteinDbRetriever.ProteomeFormat.fasta, ProteinDbRetriever.Reviewed.no,
                ProteinDbRetriever.Compress.no, ProteinDbRetriever.IncludeIsoforms.no));

            Assert.That(ex.Message, Does.Contain("UP000008595"));
            Assert.That(Directory.GetFiles(_storageDirectory), Is.Empty);
            return Task.CompletedTask;
        });
}
