using NUnit.Framework;
using Readers;
using Readers.ExternalResults.ResultFiles;
using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Security.Cryptography;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// Live canary: an mzIdentML file as PRIDE serves it, downloaded with PrideArchiveClient and read with
/// MzIdentMLResultFile straight from the .gz. It pins one small, old file rather than scanning projects.
///
/// Anything PRIDE does wrong SKIPS, via <see cref="ExternalServiceTestHelper.RunAsync"/>. That covers an
/// outage, a stalled or cut-off transfer, an empty or non-gzip body, and running past the deadline. The
/// file vanishing from the manifest or losing its HTTPS location is drift, and is reported as
/// inconclusive. The test FAILS only when intact gzip bytes do not read, which is the canary's job:
/// PRIDE hosting an mzIdentML flavour mzLib cannot read. Parser correctness is covered offline by
/// TestMzIdentMLResultFile.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideMzIdentMLLiveTests
{
    private const string Accession = "PXD000710";
    private const string FileName = "ma190_19_tandem_pproph.pep.mzid.gz";

    // as downloaded 2026-09-17: 8,563 bytes, X!Tandem + PeptideProphet, 31 SpectrumIdentificationItems
    private const string PinnedSha256 = "80bee1aa7c6637a8e780d6f7528cae0ac82e031b06f8791879d936a905cae99f";
    private const int PinnedItems = 31;

    [Test]
    public Task MzIdentMLResultFile_ReadsAPrideMzidGzAsServed() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            // a slow-but-moving transfer is bounded by nothing else, and the CI job's own timeout fails the run
            using var deadline = new CancellationTokenSource(TimeSpan.FromMinutes(3));
            using var client = new PrideArchiveClient { BodyStallTimeout = TimeSpan.FromSeconds(30) };
            string dir = Path.Combine(Path.GetTempPath(), "PrideLiveMzid", Guid.NewGuid().ToString("N"));
            try
            {
                string path;
                try
                {
                    var files = await client.GetProjectFilesAsync(Accession, cancellationToken: deadline.Token);
                    var file = files.FirstOrDefault(f => f.FileName == FileName);
                    Assume.That(file, Is.Not.Null, $"{FileName} is no longer listed for {Accession}");
                    Assume.That(file!.TryGetHttpsDownloadUrl(out _), Is.True, $"{FileName} has no HTTPS-reachable location");

                    path = await client.DownloadFileAsync(file, dir, cancellationToken: deadline.Token);
                }
                catch (OperationCanceledException) when (deadline.IsCancellationRequested)
                {
                    throw new ExternalServiceUnavailableException("PRIDE did not deliver the file within the deadline");
                }

                byte[] bytes = File.ReadAllBytes(path);
                if (bytes.Length < 2 || bytes[0] != 0x1f || bytes[1] != 0x8b)
                {
                    throw new ExternalServiceUnavailableException($"PRIDE served {bytes.Length} bytes that are not gzip");
                }

                try
                {
                    using var gzip = new GZipStream(new MemoryStream(bytes), CompressionMode.Decompress);
                    gzip.CopyTo(Stream.Null);
                }
                catch (InvalidDataException e)
                {
                    throw new ExternalServiceUnavailableException($"the transfer was not intact gzip: {e.Message}");
                }

                var mzid = new MzIdentMLResultFile(path);

                Assert.That(mzid.FileType, Is.EqualTo(SupportedFileType.MzIdentMLGz));
                Assert.That(mzid.Results, Is.Not.Empty);
                Assert.That(mzid.Results.All(r => r.OneBasedScanNumber > 0 && r.BaseSequence.Length > 0), Is.True);

                if (Convert.ToHexString(SHA256.HashData(bytes)).Equals(PinnedSha256, StringComparison.OrdinalIgnoreCase))
                {
                    Assert.That(mzid.Results, Has.Count.EqualTo(PinnedItems));
                    Assert.That(mzid.SkippedMatches, Is.Empty);
                    Assert.That(mzid.Results[0].FileNameWithoutExtension, Is.EqualTo("ma190_19"));
                }
            }
            finally
            {
                if (Directory.Exists(dir))
                {
                    Directory.Delete(dir, recursive: true);
                }
            }
        });
}
