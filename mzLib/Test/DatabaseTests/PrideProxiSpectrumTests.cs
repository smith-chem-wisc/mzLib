using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideProxiSpectrumTests
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

    private static HttpResponseMessage JsonResponse(string body, HttpStatusCode status = HttpStatusCode.OK) =>
        new(status) { Content = new StringContent(body, Encoding.UTF8, "application/json") };

    /// <summary>A PROXI spectrum object in the real wire shape (lowercase keys; cvParam attributes with no cvLabel).</summary>
    private static string SpectrumJson(string usi, double[] mzs, double[] intensities)
    {
        string Arr(double[] xs) => "[" + string.Join(",", xs.Select(x => x.ToString(CultureInfo.InvariantCulture))) + "]";
        return $$"""
        {
          "status": "READABLE",
          "usi": "{{usi}}",
          "mzs": {{Arr(mzs)}},
          "intensities": {{Arr(intensities)}},
          "attributes": [
            { "accession": "MS:1000041", "name": "charge state", "value": "2" },
            { "accession": "MS:1001910", "name": "LTQ Orbitrap Elite" }
          ]
        }
        """;
    }

    private static string ProxiArray(params string[] spectrumJson) => "[" + string.Join(",", spectrumJson) + "]";

    // ---- ToMzSpectrum (pure) -----------------------------------------------

    [Test]
    public void ToMzSpectrum_MapsPeaks()
    {
        var proxi = new PrideProxiSpectrum
        {
            Usi = "mzspec:PXD000001:file:scan:1:PEPTIDE/2",
            Mzs = new[] { 100.0, 200.0, 300.0 },
            Intensities = new[] { 10.0, 20.0, 30.0 },
        };

        MzSpectrum spectrum = proxi.ToMzSpectrum();

        Assert.That(spectrum.XArray, Is.EqualTo(new[] { 100.0, 200.0, 300.0 }));
        Assert.That(spectrum.YArray, Is.EqualTo(new[] { 10.0, 20.0, 30.0 }));
        Assert.That(spectrum.Size, Is.EqualTo(3));
    }

    [Test]
    public void ToMzSpectrum_SortsPeaksByAscendingMz_KeepingIntensityPairing()
    {
        // MzSpectrum assumes ascending m/z; ToMzSpectrum must sort defensively and carry intensities along.
        var proxi = new PrideProxiSpectrum
        {
            Mzs = new[] { 300.0, 100.0, 200.0 },
            Intensities = new[] { 30.0, 10.0, 20.0 },
        };

        MzSpectrum spectrum = proxi.ToMzSpectrum();

        Assert.That(spectrum.XArray, Is.EqualTo(new[] { 100.0, 200.0, 300.0 }));
        Assert.That(spectrum.YArray, Is.EqualTo(new[] { 10.0, 20.0, 30.0 })); // intensity 10 still paired with m/z 100
    }

    [Test]
    public void ToMzSpectrum_DoesNotReorderTheSourceDtoArrays()
    {
        // the sort works on a copy, so a caller re-reading the DTO sees its original (server) order preserved
        var proxi = new PrideProxiSpectrum
        {
            Mzs = new[] { 300.0, 100.0, 200.0 },
            Intensities = new[] { 30.0, 10.0, 20.0 },
        };

        proxi.ToMzSpectrum();

        Assert.That(proxi.Mzs, Is.EqualTo(new[] { 300.0, 100.0, 200.0 }));
        Assert.That(proxi.Intensities, Is.EqualTo(new[] { 30.0, 10.0, 20.0 }));
    }

    [Test]
    public void ToMzSpectrum_EmptyPeaks_ReturnsEmptySpectrum()
    {
        var proxi = new PrideProxiSpectrum(); // default empty arrays

        MzSpectrum spectrum = proxi.ToMzSpectrum();

        Assert.That(spectrum.Size, Is.EqualTo(0));
    }

    [Test]
    public void ToMzSpectrum_NullPeakArrays_TreatedAsEmpty()
    {
        // a deserializer could leave an array null; ToMzSpectrum coalesces to empty rather than throwing NRE
        var proxi = new PrideProxiSpectrum { Mzs = null, Intensities = null };

        MzSpectrum spectrum = proxi.ToMzSpectrum();

        Assert.That(spectrum.Size, Is.EqualTo(0));
    }

    [Test]
    public void ToMzSpectrum_LengthMismatch_ThrowsMzLibException()
    {
        var proxi = new PrideProxiSpectrum
        {
            Usi = "mzspec:PXD000001:file:scan:1:PEPTIDE/2",
            Mzs = new[] { 100.0, 200.0 },
            Intensities = new[] { 10.0 }, // one short
        };

        var ex = Assert.Throws<MzLibException>(() => proxi.ToMzSpectrum());
        Assert.That(ex.Message, Does.Contain("PEPTIDE")); // message names the offending USI
    }

    [Test]
    public void ToMzSpectrum_NullSpectrum_ThrowsArgumentNullException()
    {
        PrideProxiSpectrum proxi = null;
        Assert.That(() => proxi.ToMzSpectrum(), Throws.ArgumentNullException);
    }

    // ---- GetProxiSpectrumAsync (offline) -----------------------------------

    [Test]
    public async Task GetProxiSpectrumAsync_DeserializesSpectrum_AndAttributes()
    {
        string usi = "mzspec:PXD000561:file:scan:17555:VLHPLEGAVVIIFK/2";
        var handler = new StubHandler(_ =>
            JsonResponse(ProxiArray(SpectrumJson(usi, new[] { 110.07, 111.06 }, new[] { 39316.4, 319.7 }))));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        PrideProxiSpectrum spectrum = await client.GetProxiSpectrumAsync(usi);

        Assert.That(spectrum.Status, Is.EqualTo("READABLE"));
        Assert.That(spectrum.Usi, Is.EqualTo(usi));
        Assert.That(spectrum.Mzs, Is.EqualTo(new[] { 110.07, 111.06 }));
        Assert.That(spectrum.Intensities, Is.EqualTo(new[] { 39316.4, 319.7 }));
        Assert.That(spectrum.Attributes.Count, Is.EqualTo(2));
        // attributes map onto CvParam by accession; the charge-state cvParam carries a value, the instrument one does not
        var charge = spectrum.Attributes.FirstOrDefault(a => a.Accession == "MS:1000041");
        Assert.That(charge, Is.Not.Null);
        Assert.That(charge.Value, Is.EqualTo("2"));
        Assert.That(spectrum.Attributes.Any(a => a.Accession == "MS:1001910" && a.Name == "LTQ Orbitrap Elite"));
    }

    [Test]
    public async Task GetProxiSpectrumAsync_BuildsAbsoluteProxiUrl_WithEscapedUsi_AndResultTypeFull()
    {
        string usi = "mzspec:PXD000561:file:scan:17555:VLHPLEGAVVIIFK/2";
        var handler = new StubHandler(_ => JsonResponse(ProxiArray(SpectrumJson(usi, new[] { 100.0 }, new[] { 1.0 }))));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        await client.GetProxiSpectrumAsync(usi);

        string requested = handler.RequestedUris.Single();
        // absolute PROXI path root (NOT the v3 archive BaseAddress), USI percent-encoded, resultType=full present
        Assert.That(requested, Does.StartWith(PrideArchiveClient.DefaultProxiBaseAddress + "spectra?usi="));
        Assert.That(requested, Does.Contain(Uri.EscapeDataString(usi)));
        Assert.That(requested, Does.Contain("resultType=full"));
    }

    [Test]
    public void GetProxiSpectrumAsync_BlankUsi_ThrowsArgumentException()
    {
        using var client = new PrideArchiveClient(new HttpClient(new StubHandler(_ => JsonResponse("[]"))));
        Assert.That(async () => await client.GetProxiSpectrumAsync("   "), Throws.ArgumentException);
    }

    [TestCase(HttpStatusCode.NotFound)]   // PROXI: unknown / unreadable USI
    [TestCase(HttpStatusCode.BadRequest)] // PROXI: malformed USI
    public void GetProxiSpectrumAsync_NonSuccessStatus_ThrowsHttpRequestException(HttpStatusCode status)
    {
        var handler = new StubHandler(_ => JsonResponse("", status));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        Assert.That(async () => await client.GetProxiSpectrumAsync("mzspec:PXD:file:scan:1:PEP/2"),
            Throws.InstanceOf<HttpRequestException>());
    }

    [TestCase("[]")]     // empty array on a 200 (contract oddity)
    [TestCase("null")]   // JSON null body
    [TestCase("[null]")] // array whose single element is null
    public void GetProxiSpectrumAsync_NoSpectrumInBody_ThrowsHttpRequestException(string body)
    {
        var handler = new StubHandler(_ => JsonResponse(body));
        using var client = new PrideArchiveClient(new HttpClient(handler));
        Assert.That(async () => await client.GetProxiSpectrumAsync("mzspec:PXD:file:scan:1:PEP/2"),
            Throws.InstanceOf<HttpRequestException>());
    }

    // ---- GetSpectrumAsync (offline) ----------------------------------------

    [Test]
    public async Task GetSpectrumAsync_ReturnsMzSpectrum_WithSortedPeaks()
    {
        string usi = "mzspec:PXD000561:file:scan:17555:VLHPLEGAVVIIFK/2";
        // deliberately unsorted on the wire to prove the end-to-end call sorts
        var handler = new StubHandler(_ =>
            JsonResponse(ProxiArray(SpectrumJson(usi, new[] { 300.0, 100.0, 200.0 }, new[] { 30.0, 10.0, 20.0 }))));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        MzSpectrum spectrum = await client.GetSpectrumAsync(usi);

        Assert.That(spectrum.XArray, Is.EqualTo(new[] { 100.0, 200.0, 300.0 }));
        Assert.That(spectrum.YArray, Is.EqualTo(new[] { 10.0, 20.0, 30.0 }));
    }
}

/// <summary>
/// Live canary for the PROXI spectrum path against the real PRIDE PROXI API. Carries
/// [Category("ExternalService")] so CI runs it in the non-blocking external-service job, and
/// [Category("Pride")] to select it alone. Routed through <see cref="ExternalServiceTestHelper.RunAsync"/>
/// so a PRIDE/EBI outage SKIPS, while a genuine contract break (wrong path, unparseable response, missing
/// peaks) FAILS. Fetches a known-stable spectrum and checks it carries peaks and metadata.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideProxiSpectrumLiveTests
{
    // A long-public spectrum from PXD000561 ("A draft map of the human proteome"), the canonical PROXI example.
    private const string KnownUsi = "mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2";

    [Test]
    public Task GetProxiSpectrumAsync_LiveKnownUsi_ReturnsPeaksAndAttributes() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            var spectrum = await client.GetProxiSpectrumAsync(KnownUsi);

            Assert.That(spectrum.Usi, Is.EqualTo(KnownUsi));
            Assert.That(spectrum.Mzs.Length, Is.GreaterThan(1));
            Assert.That(spectrum.Mzs.Length, Is.EqualTo(spectrum.Intensities.Length));
            Assert.That(spectrum.Attributes, Is.Not.Empty);
        });

    [Test]
    public Task GetSpectrumAsync_LiveKnownUsi_ReturnsPopulatedMzSpectrum() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient();
            var spectrum = await client.GetSpectrumAsync(KnownUsi);

            Assert.That(spectrum.Size, Is.GreaterThan(1));
            // peaks come back m/z-ascending
            Assert.That(spectrum.XArray.First(), Is.LessThan(spectrum.XArray.Last()));
        });
}
