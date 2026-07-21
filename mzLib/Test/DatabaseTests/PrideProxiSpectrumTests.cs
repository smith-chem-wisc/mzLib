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

    private static CvParam Cv(string accession, string name, string value = "",
        string unitAccession = "", string unitName = "") =>
        new() { Accession = accession, Name = name, Value = value, UnitAccession = unitAccession, UnitName = unitName };

    /// <summary>
    /// The attribute set PRIDE really returns for a fragment spectrum, verified live against
    /// mzspec:PXD000561:Adult_Frontalcortex_bRP_Elite_85_f09:scan:17555:VLHPLEGAVVIIFK/2.
    /// </summary>
    private static List<CvParam> RealisticAttributes() => new()
    {
        Cv("MS:1001910", "LTQ Orbitrap Elite"),
        Cv("MS:1000463", "instrument", "LTQ Orbitrap Velos"),
        Cv("MS:1000463", "instrument", "LTQ Orbitrap Elite"), // repeats: one term per project instrument
        Cv("MS:1000512", "filter string", "FTMS + p NSI d Full ms2 767.97@hcd32.00 [110.00-1550.00]"),
        Cv("MS:1000828", "isolation window lower offset", "1"),
        Cv("MS:1000041", "charge state", "2"),
        Cv("MS:1000511", "ms level", "2"),
        Cv("MS:1000525", "spectrum representation", "centroid spectrum"),
        Cv("MS:1003057", "scan number", "17555"),
        Cv("MS:1000827", "isolation window target m/z", "767.973937988281"),
        Cv("MS:1000016", "scan start time", "4662.3241"),
        Cv("MS:1000829", "isolation window upper offset", "1"),
        Cv("MS:1000927", "ion injection time", "147.372"),
        Cv("MS:1000465", "scan polarity", "positive scan"),
        Cv("MS:1000744", "selected ion m/z", "767.973937988281"),
    };

    // ---- PrideProxiSpectrum is an MzSpectrum -------------------------------

    [Test]
    public void Constructor_MapsPeaksOntoTheSpectrumItself()
    {
        var proxi = new PrideProxiSpectrum(
            new[] { 100.0, 200.0, 300.0 }, new[] { 10.0, 20.0, 30.0 },
            usi: "mzspec:PXD000001:file:scan:1:PEPTIDE/2");

        Assert.That(proxi, Is.InstanceOf<MzSpectrum>()); // no conversion step needed
        Assert.That(proxi.XArray, Is.EqualTo(new[] { 100.0, 200.0, 300.0 }));
        Assert.That(proxi.YArray, Is.EqualTo(new[] { 10.0, 20.0, 30.0 }));
        Assert.That(proxi.Size, Is.EqualTo(3));
        Assert.That(proxi.Usi, Is.EqualTo("mzspec:PXD000001:file:scan:1:PEPTIDE/2"));
    }

    [Test]
    public void Constructor_SortsPeaksByAscendingMz_KeepingIntensityPairing()
    {
        // MzSpectrum assumes ascending m/z and never sorts, so the DTO must sort defensively on construction.
        var proxi = new PrideProxiSpectrum(new[] { 300.0, 100.0, 200.0 }, new[] { 30.0, 10.0, 20.0 });

        Assert.That(proxi.XArray, Is.EqualTo(new[] { 100.0, 200.0, 300.0 }));
        Assert.That(proxi.YArray, Is.EqualTo(new[] { 10.0, 20.0, 30.0 })); // intensity 10 still paired with m/z 100
    }

    [Test]
    public void Constructor_DoesNotReorderTheCallersArrays()
    {
        // the base ctor copies, so the sort cannot reach back and reorder arrays the caller still holds
        double[] mzs = { 300.0, 100.0, 200.0 };
        double[] intensities = { 30.0, 10.0, 20.0 };

        _ = new PrideProxiSpectrum(mzs, intensities);

        Assert.That(mzs, Is.EqualTo(new[] { 300.0, 100.0, 200.0 }));
        Assert.That(intensities, Is.EqualTo(new[] { 30.0, 10.0, 20.0 }));
    }

    [Test]
    public void Constructor_NoArguments_IsAnEmptySpectrumWithNoMetadata()
    {
        var proxi = new PrideProxiSpectrum();

        Assert.That(proxi.Size, Is.EqualTo(0));
        Assert.That(proxi.Usi, Is.Empty);
        Assert.That(proxi.Status, Is.Empty);
        Assert.That(proxi.Attributes, Is.Empty); // never null, so callers can enumerate unguarded
    }

    [Test]
    public void Constructor_NullPeakArrays_TreatedAsEmpty()
    {
        // a server could omit the peak arrays entirely; coalesce to empty rather than throwing NRE in the base ctor
        var proxi = new PrideProxiSpectrum(null, null, usi: "mzspec:PXD000001:file:scan:1:PEPTIDE/2");

        Assert.That(proxi.Size, Is.EqualTo(0));
    }

    [Test]
    public void Constructor_LengthMismatch_ThrowsMzLibException()
    {
        var ex = Assert.Throws<MzLibException>(() => new PrideProxiSpectrum(
            new[] { 100.0, 200.0 }, new[] { 10.0 }, // one intensity short
            usi: "mzspec:PXD000001:file:scan:1:PEPTIDE/2"));

        Assert.That(ex.Message, Does.Contain("PEPTIDE")); // message names the offending USI
    }

    // ---- ToMsDataScan (pure) -----------------------------------------------

    [Test]
    public void ToMsDataScan_ReadsEveryCvTermPrideActuallySends()
    {
        var proxi = new PrideProxiSpectrum(
            new[] { 110.07, 111.06, 1500.0 }, new[] { 39316.4, 319.7, 12.0 },
            usi: "mzspec:PXD000561:file:scan:17555:VLHPLEGAVVIIFK/2",
            status: "READABLE",
            attributes: RealisticAttributes());

        MsDataScan scan = proxi.ToMsDataScan();

        Assert.Multiple(() =>
        {
            Assert.That(scan.MassSpectrum, Is.SameAs(proxi)); // the DTO is carried through, not re-wrapped
            Assert.That(scan.OneBasedScanNumber, Is.EqualTo(17555));  // MS:1003057
            Assert.That(scan.MsnOrder, Is.EqualTo(2));                // MS:1000511
            Assert.That(scan.IsCentroid, Is.True);                    // MS:1000525 "centroid spectrum"
            Assert.That(scan.Polarity, Is.EqualTo(Polarity.Positive));// MS:1000465 "positive scan"
            Assert.That(scan.RetentionTime, Is.EqualTo(4662.3241 / 60.0).Within(1e-9)); // MS:1000016, seconds -> minutes
            Assert.That(scan.SelectedIonMZ, Is.EqualTo(767.973937988281).Within(1e-9)); // MS:1000744
            Assert.That(scan.SelectedIonChargeStateGuess, Is.EqualTo(2));               // MS:1000041
            Assert.That(scan.IsolationMz, Is.EqualTo(767.973937988281).Within(1e-9));   // MS:1000827
            Assert.That(scan.IsolationWidth, Is.EqualTo(2.0).Within(1e-9));             // MS:1000828 + MS:1000829
            Assert.That(scan.InjectionTime, Is.EqualTo(147.372).Within(1e-9));          // MS:1000927
            Assert.That(scan.ScanFilter, Does.Contain("hcd32.00"));                     // MS:1000512, verbatim
            Assert.That(scan.NativeId, Is.EqualTo("scan=17555"));
            Assert.That(scan.TotalIonCurrent, Is.EqualTo(39316.4 + 319.7 + 12.0).Within(1e-9));
            // the scan window is synthesized from the peaks, since PROXI states none
            Assert.That(scan.ScanWindowRange.Minimum, Is.EqualTo(110.07).Within(1e-9));
            Assert.That(scan.ScanWindowRange.Maximum, Is.EqualTo(1500.0).Within(1e-9));
        });
    }

    [Test]
    public void ToMsDataScan_LeavesUnstatedFieldsUnknown_RatherThanGuessing()
    {
        // PROXI names the instrument but not the analyzer, and carries no dissociation term at all
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 }, attributes: RealisticAttributes());

        MsDataScan scan = proxi.ToMsDataScan();

        Assert.That(scan.MzAnalyzer, Is.EqualTo(MZAnalyzerType.Unknown));
        Assert.That(scan.DissociationType, Is.EqualTo(DissociationType.Unknown));
        Assert.That(scan.OneBasedPrecursorScanNumber, Is.Null);
        Assert.That(scan.NoiseData, Is.Null);
    }

    [Test]
    public void ToMsDataScan_NoAttributes_FallsBackToDefaults()
    {
        var proxi = new PrideProxiSpectrum(new[] { 100.0, 200.0 }, new[] { 1.0, 2.0 });

        MsDataScan scan = proxi.ToMsDataScan();

        Assert.Multiple(() =>
        {
            Assert.That(scan.OneBasedScanNumber, Is.EqualTo(1));
            Assert.That(scan.MsnOrder, Is.EqualTo(2)); // a USI addresses a fragment spectrum
            Assert.That(scan.IsCentroid, Is.False);
            Assert.That(scan.Polarity, Is.EqualTo(Polarity.Unknown));
            Assert.That(scan.RetentionTime, Is.EqualTo(0));
            Assert.That(scan.SelectedIonMZ, Is.Null);
            Assert.That(scan.SelectedIonChargeStateGuess, Is.Null);
            Assert.That(scan.IsolationMz, Is.Null);
            Assert.That(scan.IsolationWidth, Is.Null);
            Assert.That(scan.InjectionTime, Is.Null);
            Assert.That(scan.ScanFilter, Is.Null);
        });
    }

    [Test]
    public void ToMsDataScan_EmptyPeaks_HasNoScanWindow()
    {
        MsDataScan scan = new PrideProxiSpectrum().ToMsDataScan();

        Assert.That(scan.ScanWindowRange, Is.Null); // nothing to bound it with
        Assert.That(scan.TotalIonCurrent, Is.EqualTo(0));
    }

    [Test]
    public void ToMsDataScan_ReadsPolarityAndCentroidFromBareTerms()
    {
        // a PROXI server other than PRIDE may send the PSI terms directly instead of as a value
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 }, attributes: new List<CvParam>
        {
            Cv("MS:1000129", "negative scan"),
            Cv("MS:1000127", "centroid spectrum"),
        });

        MsDataScan scan = proxi.ToMsDataScan();

        Assert.That(scan.Polarity, Is.EqualTo(Polarity.Negative));
        Assert.That(scan.IsCentroid, Is.True);
    }

    [Test]
    public void ToMsDataScan_UnrecognizedPolarityValue_IsUnknown()
    {
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 },
            attributes: new List<CvParam> { Cv("MS:1000465", "scan polarity", "sideways scan") });

        Assert.That(proxi.ToMsDataScan().Polarity, Is.EqualTo(Polarity.Unknown));
    }

    [TestCase("UO:0000031", "minute", 12.5)]  // stated in minutes -> taken as-is
    [TestCase("UO:0000010", "second", 750.0 / 60.0)] // stated in seconds -> converted
    [TestCase("", "", 750.0 / 60.0)]          // no unit -> seconds, which is what PRIDE serves
    public void ToMsDataScan_NormalizesRetentionTimeToMinutes(string unitAccession, string unitName, double expected)
    {
        double value = unitName == "minute" ? 12.5 : 750.0;
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 }, attributes: new List<CvParam>
        {
            Cv("MS:1000016", "scan start time", value.ToString(CultureInfo.InvariantCulture), unitAccession, unitName),
        });

        Assert.That(proxi.ToMsDataScan().RetentionTime, Is.EqualTo(expected).Within(1e-9));
    }

    [Test]
    public void ToMsDataScan_IsolationMzFallsBackToPrecursorMz_WhenNoWindowStated()
    {
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 },
            attributes: new List<CvParam> { Cv("MS:1000744", "selected ion m/z", "500.25") });

        MsDataScan scan = proxi.ToMsDataScan();

        Assert.That(scan.IsolationMz, Is.EqualTo(500.25).Within(1e-9));
        Assert.That(scan.IsolationWidth, Is.Null); // no offsets stated, so no width is invented
    }

    [Test]
    public void ToMsDataScan_OneIsolationOffsetOnly_UsesIt()
    {
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 }, attributes: new List<CvParam>
        {
            Cv("MS:1000827", "isolation window target m/z", "500.0"),
            Cv("MS:1000829", "isolation window upper offset", "1.5"),
        });

        Assert.That(proxi.ToMsDataScan().IsolationWidth, Is.EqualTo(1.5).Within(1e-9));
    }

    [Test]
    public void ToMsDataScan_NonNumericCvValue_IsIgnoredRatherThanThrowing()
    {
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 }, attributes: new List<CvParam>
        {
            Cv("MS:1000041", "charge state", "not a number"),
            Cv("MS:1000511", "ms level", ""),
        });

        MsDataScan scan = proxi.ToMsDataScan();

        Assert.That(scan.SelectedIonChargeStateGuess, Is.Null);
        Assert.That(scan.MsnOrder, Is.EqualTo(2)); // falls back to the default
    }

    [Test]
    public void ToMsDataScan_ReadsTheFirstOccurrenceOfARepeatedTerm()
    {
        // PRIDE sends one "instrument" term per project instrument; a repeated term must not throw
        var proxi = new PrideProxiSpectrum(new[] { 100.0 }, new[] { 1.0 }, attributes: new List<CvParam>
        {
            Cv("MS:1003057", "scan number", "111"),
            Cv("MS:1003057", "scan number", "222"),
        });

        Assert.That(proxi.ToMsDataScan().OneBasedScanNumber, Is.EqualTo(111));
    }

    [Test]
    public void ToMsDataScan_NullSpectrum_ThrowsArgumentNullException()
    {
        PrideProxiSpectrum proxi = null;
        Assert.That(() => proxi.ToMsDataScan(), Throws.ArgumentNullException);
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
        Assert.That(spectrum.XArray, Is.EqualTo(new[] { 110.07, 111.06 }));
        Assert.That(spectrum.YArray, Is.EqualTo(new[] { 39316.4, 319.7 }));
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

    [Test]
    public void GetProxiSpectrumAsync_RaggedPeakArraysOnTheWire_ThrowsMzLibException()
    {
        // the parallel-arrays guard runs inside the constructor Newtonsoft deserializes through; it surfaces
        // unwrapped rather than as a JsonSerializationException, so the documented contract still holds
        var handler = new StubHandler(_ => JsonResponse(
            """[{"usi":"mzspec:PXD000001:file:scan:1:PEPTIDE/2","mzs":[100.0,200.0],"intensities":[1.0]}]"""));
        using var client = new PrideArchiveClient(new HttpClient(handler));

        var ex = Assert.ThrowsAsync<MzLibException>(
            async () => await client.GetProxiSpectrumAsync("mzspec:PXD000001:file:scan:1:PEPTIDE/2"));
        Assert.That(ex.Message, Does.Contain("PEPTIDE"));
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
            Assert.That(spectrum.Size, Is.GreaterThan(1));
            Assert.That(spectrum.XArray.Length, Is.EqualTo(spectrum.YArray.Length));
            Assert.That(spectrum.Attributes, Is.Not.Empty);

            // the CV attributes read into a real MsDataScan
            MsDataScan scan = spectrum.ToMsDataScan();
            Assert.That(scan.MsnOrder, Is.EqualTo(2));
            Assert.That(scan.OneBasedScanNumber, Is.EqualTo(17555)); // the scan number in the USI
            Assert.That(scan.SelectedIonChargeStateGuess, Is.EqualTo(2)); // the charge in the USI
            Assert.That(scan.SelectedIonMZ, Is.Not.Null);
            Assert.That(scan.Polarity, Is.EqualTo(Polarity.Positive));
            Assert.That(scan.RetentionTime, Is.GreaterThan(0));
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
