using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Collections.Generic;
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
public class PrideProjectTests
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

    private static HttpResponseMessage JsonResponse(string body, HttpStatusCode status = HttpStatusCode.OK)
        => new(status) { Content = new StringContent(body, Encoding.UTF8, "application/json") };

    /// <summary>
    /// A PRIDE project metadata object in the real v3 wire shape, captured live from
    /// <c>projects/PXD012345</c> on 2026-07-21 and trimmed to one entry per collection. Note the
    /// bare calendar dates, the "@type" discriminators, and that sampleAttributes nests a CvParam
    /// key against a list of CvParam values.
    /// </summary>
    private const string ProjectJson =
        """
        {
          "accession": "PXD012345",
          "title": "Metaproteomics benchmark dataset",
          "projectDescription": "A benchmark of four organisms.",
          "sampleProcessingProtocol": "Cells were lysed and digested with trypsin.",
          "dataProcessingProtocol": "Searched with MaxQuant against a combined database.",
          "doi": "10.6019/PXD012345",
          "submissionType": "COMPLETE",
          "license": "CC BY 4.0",
          "submissionDate": "2019-01-15",
          "publicationDate": "2025-07-13",
          "projectTags": ["Technical", "Metaproteomics"],
          "keywords": ["benchmark"],
          "countries": ["Italy"],
          "submitters": [
            {
              "title": "Professor",
              "firstName": "Anna Laura",
              "lastName": "Capriotti",
              "identifier": "11200124",
              "affiliation": "Sapienza University of Rome",
              "email": "annalaura.capriotti@uniroma1.it",
              "country": "Italy",
              "orcid": "",
              "name": "Anna Laura Capriotti",
              "id": "11200124"
            }
          ],
          "labPIs": [
            { "title": "Dr", "firstName": "Aldo", "lastName": "Lagana", "affiliation": "Sapienza", "email": "", "country": "Italy", "orcid": "0000-0002-1234-5678", "name": "Aldo Lagana", "id": "99" }
          ],
          "references": [
            { "referenceLine": "Macedo-Silva C et al. Cell Death Discov. 2025 11(1):306", "pubmedID": 40610411, "doi": "10.1038/s41420-025-02597-4" }
          ],
          "instruments": [ { "@type": "CvParam", "cvLabel": "MS", "accession": "MS:1000449", "name": "LTQ Orbitrap" } ],
          "softwares": [ { "@type": "CvParam", "cvLabel": "MS", "accession": "MS:1001583", "name": "MaxQuant" } ],
          "experimentTypes": [ { "@type": "CvParam", "cvLabel": "MS", "accession": "MS:1000550", "name": "Shotgun proteomics" } ],
          "quantificationMethods": [],
          "organisms": [ { "@type": "CvParam", "cvLabel": "NEWT", "accession": "NEWT:4530", "name": "Oryza sativa (rice)" } ],
          "organismParts": [],
          "diseases": [],
          "identifiedPTMStrings": [ { "@type": "CvParam", "cvLabel": "PRIDE", "accession": "PRIDE:0000398", "name": "No PTMs are included in the dataset" } ],
          "additionalAttributes": [],
          "sampleAttributes": [
            {
              "@type": "Tuple",
              "key": { "cvLabel": "EFO", "accession": "OBI:0100026", "name": "organism" },
              "value": [
                { "cvLabel": "NEWT", "accession": "NEWT:5476", "name": "Candida albicans (Yeast)" },
                { "cvLabel": "NEWT", "accession": "NEWT:562", "name": "Escherichia coli" }
              ]
            }
          ],
          "totalFileDownloads": 19,
          "botCount": 3,
          "hubCount": 1,
          "organicCount": 15
        }
        """;

    private static PrideArchiveClient ClientReturning(string body, HttpStatusCode status = HttpStatusCode.OK,
        StubHandler handler = null)
    {
        handler ??= new StubHandler(_ => JsonResponse(body, status));
        return new PrideArchiveClient(new HttpClient(handler));
    }

    // ---- deserialization ----------------------------------------------------

    [Test]
    public async Task GetProjectAsync_RealWireShape_DeserializesScalarFields()
    {
        using var client = ClientReturning(ProjectJson);

        PrideProject project = await client.GetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(project.Accession, Is.EqualTo("PXD012345"));
            Assert.That(project.Title, Is.EqualTo("Metaproteomics benchmark dataset"));
            Assert.That(project.ProjectDescription, Is.EqualTo("A benchmark of four organisms."));
            Assert.That(project.SampleProcessingProtocol, Does.StartWith("Cells were lysed"));
            Assert.That(project.DataProcessingProtocol, Does.StartWith("Searched with MaxQuant"));
            Assert.That(project.Doi, Is.EqualTo("10.6019/PXD012345"));
            Assert.That(project.SubmissionType, Is.EqualTo("COMPLETE"));
            Assert.That(project.License, Is.EqualTo("CC BY 4.0"));
            Assert.That(project.TotalFileDownloads, Is.EqualTo(19));
            Assert.That(project.BotCount, Is.EqualTo(3));
            Assert.That(project.HubCount, Is.EqualTo(1));
            Assert.That(project.OrganicCount, Is.EqualTo(15));
        });
    }

    /// <summary>
    /// The dates arrive as bare calendar dates with no time and no UTC offset, which is why they are
    /// modeled as DateTime rather than the DateTimeOffset used on PrideArchiveFile.
    /// </summary>
    /// <remarks>
    /// DateTime equality ignores Kind, so comparing values alone would still pass if the deserializer
    /// attached a machine-local offset — the exact drift this test exists to catch. Kind is therefore
    /// asserted explicitly: Unspecified means no zone was invented, and the time component must be
    /// midnight because the wire carries no time at all.
    /// </remarks>
    [Test]
    public async Task GetProjectAsync_BareCalendarDates_ParseWithoutTimeZoneDrift()
    {
        using var client = ClientReturning(ProjectJson);

        PrideProject project = await client.GetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(project.SubmissionDate, Is.EqualTo(new DateTime(2019, 1, 15)));
            Assert.That(project.PublicationDate, Is.EqualTo(new DateTime(2025, 7, 13)));
            Assert.That(project.SubmissionDate.Kind, Is.EqualTo(DateTimeKind.Unspecified));
            Assert.That(project.PublicationDate.Kind, Is.EqualTo(DateTimeKind.Unspecified));
            Assert.That(project.SubmissionDate.TimeOfDay, Is.EqualTo(TimeSpan.Zero));
            Assert.That(project.PublicationDate.TimeOfDay, Is.EqualTo(TimeSpan.Zero));
        });
    }

    [Test]
    public async Task GetProjectAsync_StringArrays_Deserialize()
    {
        using var client = ClientReturning(ProjectJson);

        PrideProject project = await client.GetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(project.ProjectTags, Is.EqualTo(new[] { "Technical", "Metaproteomics" }));
            Assert.That(project.Keywords, Is.EqualTo(new[] { "benchmark" }));
            Assert.That(project.Countries, Is.EqualTo(new[] { "Italy" }));
        });
    }

    [Test]
    public async Task GetProjectAsync_CvParamGroups_DeserializeAndKeyOnAccession()
    {
        using var client = ClientReturning(ProjectJson);

        PrideProject project = await client.GetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(project.Instruments.Single().Accession, Is.EqualTo("MS:1000449"));
            Assert.That(project.Instruments.Single().Name, Is.EqualTo("LTQ Orbitrap"));
            Assert.That(project.Instruments.Single().CvLabel, Is.EqualTo("MS"));
            Assert.That(project.Softwares.Single().Accession, Is.EqualTo("MS:1001583"));
            Assert.That(project.ExperimentTypes.Single().Accession, Is.EqualTo("MS:1000550"));
            Assert.That(project.Organisms.Single().Accession, Is.EqualTo("NEWT:4530"));
            Assert.That(project.IdentifiedPTMStrings.Single().Accession, Is.EqualTo("PRIDE:0000398"));
        });
    }

    /// <summary>Empty CV collections must arrive as empty lists, never null — the class promises "never null".</summary>
    [Test]
    public async Task GetProjectAsync_EmptyCvCollections_AreEmptyNotNull()
    {
        using var client = ClientReturning(ProjectJson);

        PrideProject project = await client.GetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(project.QuantificationMethods, Is.Empty);
            Assert.That(project.OrganismParts, Is.Empty);
            Assert.That(project.Diseases, Is.Empty);
            Assert.That(project.AdditionalAttributes, Is.Empty);
        });
    }

    [Test]
    public async Task GetProjectAsync_Submitters_DeserializeEveryContactField()
    {
        using var client = ClientReturning(ProjectJson);

        PrideContact submitter = (await client.GetProjectAsync("PXD012345")).Submitters.Single();

        Assert.Multiple(() =>
        {
            Assert.That(submitter.Title, Is.EqualTo("Professor"));
            Assert.That(submitter.FirstName, Is.EqualTo("Anna Laura"));
            Assert.That(submitter.LastName, Is.EqualTo("Capriotti"));
            Assert.That(submitter.Name, Is.EqualTo("Anna Laura Capriotti"));
            Assert.That(submitter.Affiliation, Is.EqualTo("Sapienza University of Rome"));
            Assert.That(submitter.Email, Is.EqualTo("annalaura.capriotti@uniroma1.it"));
            Assert.That(submitter.Country, Is.EqualTo("Italy"));
            Assert.That(submitter.Orcid, Is.Empty);
            Assert.That(submitter.Id, Is.EqualTo("11200124"));
        });
    }

    [Test]
    public async Task GetProjectAsync_LabPIs_Deserialize()
    {
        using var client = ClientReturning(ProjectJson);

        PrideContact pi = (await client.GetProjectAsync("PXD012345")).LabPIs.Single();

        Assert.Multiple(() =>
        {
            Assert.That(pi.Name, Is.EqualTo("Aldo Lagana"));
            Assert.That(pi.Orcid, Is.EqualTo("0000-0002-1234-5678"));
            Assert.That(pi.Email, Is.Empty);
        });
    }

    [Test]
    public async Task GetProjectAsync_References_DeserializeIncludingNumericPubmedId()
    {
        using var client = ClientReturning(ProjectJson);

        PrideReference reference = (await client.GetProjectAsync("PXD012345")).References.Single();

        Assert.Multiple(() =>
        {
            Assert.That(reference.PubmedId, Is.EqualTo(40610411));
            Assert.That(reference.Doi, Is.EqualTo("10.1038/s41420-025-02597-4"));
            Assert.That(reference.ReferenceLine, Does.Contain("Cell Death Discov"));
        });
    }

    /// <summary>
    /// sampleAttributes is the only nested shape PRIDE sends: a single CvParam key against a LIST of
    /// CvParam values. The "@type": "Tuple" discriminator must be ignored, not bound.
    /// </summary>
    [Test]
    public async Task GetProjectAsync_SampleAttributes_DeserializeKeyAndValueList()
    {
        using var client = ClientReturning(ProjectJson);

        PrideSampleAttribute attribute = (await client.GetProjectAsync("PXD012345")).SampleAttributes.Single();

        Assert.Multiple(() =>
        {
            Assert.That(attribute.Key.Accession, Is.EqualTo("OBI:0100026"));
            Assert.That(attribute.Key.Name, Is.EqualTo("organism"));
            Assert.That(attribute.Value, Has.Count.EqualTo(2));
            Assert.That(attribute.Value.Select(v => v.Accession), Is.EqualTo(new[] { "NEWT:5476", "NEWT:562" }));
        });
    }

    /// <summary>
    /// PRIDE omits fields for sparsely-annotated projects. Absent members must fall back to the DTO's
    /// non-null defaults rather than becoming null.
    /// </summary>
    [Test]
    public async Task GetProjectAsync_MinimalPayload_LeavesDefaultsIntact()
    {
        using var client = ClientReturning("""{ "accession": "PXD000001" }""");

        PrideProject project = await client.GetProjectAsync("PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(project.Accession, Is.EqualTo("PXD000001"));
            Assert.That(project.Title, Is.Empty);
            Assert.That(project.Doi, Is.Empty);
            Assert.That(project.Instruments, Is.Empty);
            Assert.That(project.Submitters, Is.Empty);
            Assert.That(project.References, Is.Empty);
            Assert.That(project.SampleAttributes, Is.Empty);
            Assert.That(project.ProjectTags, Is.Empty);
        });
    }

    /// <summary>An explicit JSON null must not clobber a non-null DTO default (NullValueHandling.Ignore).</summary>
    [Test]
    public async Task GetProjectAsync_ExplicitJsonNulls_DoNotClobberDefaults()
    {
        using var client = ClientReturning("""
            { "accession": "PXD000001", "title": null, "instruments": null, "submitters": null, "projectTags": null }
            """);

        PrideProject project = await client.GetProjectAsync("PXD000001");

        Assert.Multiple(() =>
        {
            Assert.That(project.Title, Is.Empty);
            Assert.That(project.Instruments, Is.Not.Null.And.Empty);
            Assert.That(project.Submitters, Is.Not.Null.And.Empty);
            Assert.That(project.ProjectTags, Is.Not.Null.And.Empty);
        });
    }

    // ---- request construction -----------------------------------------------

    /// <summary>
    /// This endpoint shares the v3 base address, so the request must resolve as a RELATIVE URI —
    /// unlike PROXI, which sits under a different path root and needs an absolute one.
    /// </summary>
    [Test]
    public async Task GetProjectAsync_BuildsRelativeV3Uri()
    {
        var handler = new StubHandler(_ => JsonResponse(ProjectJson));
        using var client = ClientReturning(null, handler: handler);

        await client.GetProjectAsync("PXD012345");

        Assert.That(handler.RequestedUris.Single(),
            Is.EqualTo("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD012345"));
    }

    /// <summary>
    /// The accession is caller-supplied, so it is escaped before being interpolated into the path.
    /// A "#" is the assertable case: unescaped it would truncate the request into a fragment and
    /// silently fetch the wrong project. (Uri.ToString() renders %20 and %2F back in decoded form,
    /// so those cannot be asserted through the recorded URI even though they are escaped too.)
    /// </summary>
    [Test]
    public async Task GetProjectAsync_EscapesAccessionIntoPath()
    {
        var handler = new StubHandler(_ => JsonResponse(ProjectJson));
        using var client = ClientReturning(null, handler: handler);

        await client.GetProjectAsync("PXD#frag");

        Assert.That(handler.RequestedUris.Single(),
            Is.EqualTo("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD%23frag"));
    }

    // ---- error contract -----------------------------------------------------
    // Kept offline deliberately: the live fixture routes through ExternalServiceTestHelper, which
    // converts HttpRequestException into Assert.Ignore — so a LIVE negative test would be silently
    // skipped rather than passed, and would prove nothing.

    /// <summary>
    /// InstanceOf, not Throws.ArgumentException: the latter is an exact-type constraint, which would pin
    /// the null case away from the BCL's ArgumentNullException idiom should the guard ever adopt it.
    /// Asserting no request was recorded proves the guard runs BEFORE the network call, not after.
    /// </summary>
    [Test]
    public void GetProjectAsync_BlankAccession_ThrowsArgumentExceptionWithoutCallingTheApi(
        [Values(null, "", "   ")] string accession)
    {
        var handler = new StubHandler(_ => JsonResponse(ProjectJson));
        using var client = ClientReturning(null, handler: handler);

        Assert.That(async () => await client.GetProjectAsync(accession), Throws.InstanceOf<ArgumentException>());
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    [Test]
    public void TryGetProjectAsync_BlankAccession_ThrowsArgumentExceptionWithoutCallingTheApi(
        [Values(null, "", "   ")] string accession)
    {
        var handler = new StubHandler(_ => JsonResponse(ProjectJson));
        using var client = ClientReturning(null, handler: handler);

        Assert.That(async () => await client.TryGetProjectAsync(accession), Throws.InstanceOf<ArgumentException>());
        Assert.That(handler.RequestedUris, Is.Empty);
    }

    /// <summary>An unknown accession 404s on this endpoint — unlike the manifest, which answers 200 [].</summary>
    [Test]
    public async Task TryGetProjectAsync_UnknownAccession_ReturnsFalseAndNull()
    {
        using var client = ClientReturning("", HttpStatusCode.NotFound);

        (bool found, PrideProject project) = await client.TryGetProjectAsync("PXD999999999");

        Assert.Multiple(() =>
        {
            Assert.That(found, Is.False);
            Assert.That(project, Is.Null);
        });
    }

    [Test]
    public async Task TryGetProjectAsync_KnownAccession_ReturnsTrueAndProject()
    {
        using var client = ClientReturning(ProjectJson);

        (bool found, PrideProject project) = await client.TryGetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(found, Is.True);
            Assert.That(project, Is.Not.Null);
            Assert.That(project.Accession, Is.EqualTo("PXD012345"));
        });
    }

    /// <summary>
    /// "No such project" must NOT be an HttpRequestException. That type means "the service is
    /// unavailable" and ExternalServiceTestHelper converts it to a skipped test, so a withdrawn
    /// accession would silently pass instead of failing. It is an MzLibException, and this test pins
    /// that it is not an HttpRequestException so the distinction cannot regress.
    /// </summary>
    [Test]
    public void GetProjectAsync_UnknownAccession_ThrowsMzLibExceptionNotHttpRequestException()
    {
        using var client = ClientReturning("", HttpStatusCode.NotFound);

        MzLibException exception = Assert.ThrowsAsync<MzLibException>(
            async () => await client.GetProjectAsync("PXD999999999"));

        Assert.Multiple(() =>
        {
            Assert.That(exception.Message, Does.Contain("PXD999999999"));
            Assert.That(exception, Is.Not.InstanceOf<HttpRequestException>());
        });
    }

    /// <summary>
    /// The whole point of the Try contract: only a 404 is an absence. A service outage must still
    /// throw, so that it can never be misreported to the caller as "no such project".
    /// </summary>
    [Test]
    public void TryGetProjectAsync_ServiceUnavailable_ThrowsRatherThanReportingNotFound(
        [Values(HttpStatusCode.InternalServerError, HttpStatusCode.ServiceUnavailable,
                HttpStatusCode.RequestTimeout, HttpStatusCode.TooManyRequests)] HttpStatusCode status)
    {
        using var client = ClientReturning("upstream failure", status);

        Assert.That(async () => await client.TryGetProjectAsync("PXD012345"),
            Throws.InstanceOf<HttpRequestException>());
    }

    [Test]
    public void GetProjectAsync_NonSuccessStatus_ThrowsHttpRequestException()
    {
        using var client = ClientReturning("server error", HttpStatusCode.InternalServerError);

        Assert.That(async () => await client.GetProjectAsync("PXD012345"),
            Throws.InstanceOf<HttpRequestException>());
    }

    /// <summary>
    /// A 200 carrying no project is a broken contract, not an absence — so it throws, and as an
    /// MzLibException rather than an HttpRequestException so a live run fails instead of skipping.
    /// Both the literal "null" payload and a genuinely empty body deserialize to null and must behave
    /// identically; the empty body is the case a real proxy or truncated response produces.
    /// </summary>
    [Test]
    public void TryGetProjectAsync_NoProjectOn200_ThrowsMzLibException(
        [Values("null", "", "   ")] string body)
    {
        using var client = ClientReturning(body);

        Assert.That(async () => await client.TryGetProjectAsync("PXD012345"),
            Throws.InstanceOf<MzLibException>());
    }

    /// <summary>
    /// An empty JSON object deserializes to a fully-defaulted PrideProject, which would otherwise be
    /// returned as a successful result carrying nothing. A real project always has its accession.
    /// </summary>
    [Test]
    public void TryGetProjectAsync_EmptyJsonObject_ThrowsMzLibException()
    {
        using var client = ClientReturning("{}");

        Assert.That(async () => await client.TryGetProjectAsync("PXD012345"),
            Throws.InstanceOf<MzLibException>().With.Message.Contains("no accession"));
    }

    /// <summary>
    /// NullValueHandling.Ignore suppresses null values, not null ELEMENTS, so a null inside a CV array
    /// would otherwise survive into a collection documented as never-null and NRE at first dereference.
    /// </summary>
    [Test]
    public async Task TryGetProjectAsync_NullElementsInArrays_AreDropped()
    {
        using var client = ClientReturning("""
            {
              "accession": "PXD012345",
              "instruments": [ null, { "accession": "MS:1000449", "name": "LTQ Orbitrap" }, null ],
              "projectTags": [ null, "Technical" ],
              "submitters": [ null ],
              "sampleAttributes": [
                null,
                { "key": { "accession": "OBI:0100026" }, "value": [ null, { "accession": "NEWT:562" } ] }
              ]
            }
            """);

        (_, PrideProject project) = await client.TryGetProjectAsync("PXD012345");

        Assert.Multiple(() =>
        {
            Assert.That(project.Instruments, Has.Count.EqualTo(1));
            Assert.That(project.Instruments.Single().Accession, Is.EqualTo("MS:1000449"));
            Assert.That(project.ProjectTags, Is.EqualTo(new[] { "Technical" }));
            Assert.That(project.Submitters, Is.Empty);
            Assert.That(project.SampleAttributes, Has.Count.EqualTo(1));
            Assert.That(project.SampleAttributes.Single().Value.Single().Accession, Is.EqualTo("NEWT:562"));
        });
    }

    [Test]
    public void GetProjectAsync_Cancelled_ThrowsOperationCanceledException()
    {
        using var client = ClientReturning(ProjectJson);
        using var cts = new CancellationTokenSource();
        cts.Cancel();

        Assert.That(async () => await client.GetProjectAsync("PXD012345", cts.Token),
            Throws.InstanceOf<OperationCanceledException>());
    }
}

/// <summary>
/// Live canary against the real PRIDE Archive API. Carries [Category("ExternalService")] so CI runs
/// it in the dedicated, non-blocking external-service job (see .github/workflows/dotnet.yml) rather
/// than the required unit-test run; [Category("Pride")] allows selecting it on its own. The call is
/// routed through <see cref="ExternalServiceTestHelper.RunAsync"/> so a PRIDE outage (timeout / 5xx /
/// unreachable) SKIPS the test, while a genuine contract break (wrong URL, unparseable response,
/// missing value) still FAILS. Run on demand — or in the external-service job — to detect API drift.
/// </summary>
/// <remarks>
/// What the live unknown-accession canary below can and cannot catch, stated precisely because the
/// distinction is easy to get wrong: if PRIDE ever answered a nonexistent accession with 200 and a
/// project, TryGetProjectAsync would report Found = true and the test FAILS, which is the drift worth
/// catching. If PRIDE instead answered 5xx, that is indistinguishable from an outage and RunAsync
/// skips — so the test cannot prove 404 is still 404. The exhaustive error-contract coverage
/// therefore lives offline in <see cref="PrideProjectTests"/>, where the status code is dictated by a
/// stub rather than by EBI's mood.
/// </remarks>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideProjectLiveTests
{
    [Test]
    public Task GetProjectAsync_LivePxd012345_ReturnsPopulatedMetadata() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            PrideProject project = await client.GetProjectAsync("PXD012345");

            Assert.Multiple(() =>
            {
                Assert.That(project.Accession, Is.EqualTo("PXD012345"));
                Assert.That(project.Title, Is.Not.Empty);
                Assert.That(project.SubmissionDate, Is.Not.EqualTo(default(DateTime)));
                Assert.That(project.Instruments, Is.Not.Empty);
                // Constraint form, not Assert.That(bool): a collapsed LINQ predicate reports only
                // "Expected: True, But was: False", which cannot be triaged from a CI log.
                Assert.That(project.Instruments, Is.All.Matches<CvParam>(i => !string.IsNullOrEmpty(i.Accession)));
                Assert.That(project.Submitters, Is.Not.Empty);
            });
        });

    [Test]
    public Task TryGetProjectAsync_LiveUnknownAccession_ReportsNotFound() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            (bool found, PrideProject project) = await client.TryGetProjectAsync("PXD999999999");

            Assert.Multiple(() =>
            {
                Assert.That(found, Is.False);
                Assert.That(project, Is.Null);
            });
        });
}
