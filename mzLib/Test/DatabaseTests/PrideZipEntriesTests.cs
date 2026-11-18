using MzLibUtil;
using NUnit.Framework;
using System;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.IO.Hashing;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// S4 + SDRF-P14: a zip's entries listed, and picked entries fetched, by HTTP range requests on its central directory
/// and its entries' bytes, without downloading the zip. Offline: every answer comes from <see cref="RangeServer"/>, which
/// answers ranges as PRIDE's FTP host was measured to on 2026-11-18, and backoffs are recorded instead of slept through.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class PrideZipEntriesTests
{
    // ---- test doubles -------------------------------------------------------

    /// <summary>
    /// A file on a server that answers byte ranges as PRIDE's FTP host does (measured 2026-11-18): a suffix or bounded
    /// <c>Range</c> gets a <c>206</c> with <c>Content-Range</c> and the file's <c>ETag</c>; a suffix longer than the file
    /// gets all of it as a <c>206</c>; a range past the end gets a <c>416</c>; an <c>If-Range</c> that does not match gets
    /// a <c>200</c> with the whole file. The file may be sparse: bytes outside its segments read as zero, so a zip can
    /// sit at an offset beyond 4 GiB without 4 GiB of memory.
    /// </summary>
    private sealed class RangeServer : HttpMessageHandler
    {
        private readonly List<(long Offset, byte[] Bytes)> _segments = new();

        public RangeServer(byte[] file) : this(file.Length, (0, file)) { }

        public RangeServer(long total, params (long Offset, byte[] Bytes)[] segments)
        {
            Total = total;
            _segments.AddRange(segments);
        }

        public long Total { get; }

        /// <summary>The file's ETag, quoted; null for a server that sends no validator.</summary>
        public string ETag { get; set; } = "\"zip-v1\"";

        /// <summary>Every request's <c>Range</c> and <c>If-Range</c>, in order.</summary>
        public List<(string Range, string IfRange)> Requests { get; } = new();

        /// <summary>Answers the Nth request (1-based) in place of the file, when it returns non-null.</summary>
        public Func<HttpRequestMessage, int, HttpResponseMessage> Before { get; set; }

        /// <summary>Changes the bytes of the Nth request's answer (request number, first offset, bytes) before they are sent.</summary>
        public Action<int, long, byte[]> Damage { get; set; }

        public byte[] Read(long from, long to)
        {
            var bytes = new byte[to - from + 1];
            foreach ((long offset, byte[] segment) in _segments)
            {
                long start = Math.Max(from, offset);
                long end = Math.Min(to, offset + segment.Length - 1);
                if (start <= end)
                    Array.Copy(segment, start - offset, bytes, start - from, end - start + 1);
            }
            return bytes;
        }

        protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            Requests.Add((request.Headers.Range?.ToString(), request.Headers.IfRange?.ToString()));
            int n = Requests.Count;
            if (Before?.Invoke(request, n) is HttpResponseMessage scripted)
                return Task.FromResult(scripted);

            RangeItemHeaderValue item = request.Headers.Range?.Ranges.Single();
            bool stale = request.Headers.IfRange is { } ifRange && ifRange.EntityTag?.Tag != ETag;
            if (item == null || stale)
            {
                if (Total > 10_000_000)
                    throw new InvalidOperationException("The test server would have sent the whole of a large file.");
                return Task.FromResult(Answer(HttpStatusCode.OK, Read(0, Total - 1), null));
            }

            long from = item.From ?? Math.Max(0, Total - item.To.Value);
            long to = item.From.HasValue ? Math.Min(item.To ?? Total - 1, Total - 1) : Total - 1;
            if (from >= Total)
                return Task.FromResult(new HttpResponseMessage(HttpStatusCode.RequestedRangeNotSatisfiable) { Content = new ByteArrayContent(Array.Empty<byte>()) });

            byte[] bytes = Read(from, to);
            Damage?.Invoke(n, from, bytes);
            return Task.FromResult(Answer(HttpStatusCode.PartialContent, bytes, new ContentRangeHeaderValue(from, to, Total)));
        }

        private HttpResponseMessage Answer(HttpStatusCode status, byte[] bytes, ContentRangeHeaderValue range)
        {
            var content = new ByteArrayContent(bytes);
            content.Headers.ContentRange = range;
            var response = new HttpResponseMessage(status) { Content = content };
            if (ETag != null)
                response.Headers.ETag = EntityTagHeaderValue.Parse(ETag);
            return response;
        }
    }

    /// <summary>A body that fails the test if anyone reads it: proof that a response was dropped unread.</summary>
    private sealed class UnreadableStream : Stream
    {
        public override int Read(byte[] buffer, int offset, int count) => throw new AssertionException("This body must never be read.");
        public override Task<int> ReadAsync(byte[] buffer, int offset, int count, CancellationToken cancellationToken) => throw new AssertionException("This body must never be read.");
        public override ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default) => throw new AssertionException("This body must never be read.");
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

    /// <summary>A body that delivers nothing until cancelled, as a stalled connection does.</summary>
    private sealed class StallingStream : Stream
    {
        public override int Read(byte[] buffer, int offset, int count) => throw new NotSupportedException();
        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            await Task.Delay(Timeout.Infinite, cancellationToken);
            return 0;
        }
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

    /// <summary>Hides a stream's ability to seek, so <see cref="ZipArchive"/> writes data descriptors (flag bit 3) as it would to a network stream.</summary>
    private sealed class WriteOnlyStream(Stream inner) : Stream
    {
        public override void Write(byte[] buffer, int offset, int count) => inner.Write(buffer, offset, count);
        public override bool CanRead => false;
        public override bool CanSeek => false;
        public override bool CanWrite => true;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }
        public override void Flush() => inner.Flush();
        public override int Read(byte[] buffer, int offset, int count) => throw new NotSupportedException();
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
    }

    // ---- zips -----------------------------------------------------------------

    private const string ZipUrl = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD000001/combined.zip";

    /// <summary>MaxQuant-like text: compressible, as summary.txt is.</summary>
    private static readonly byte[] Summary = Encoding.UTF8.GetBytes(string.Join("\n",
        Enumerable.Range(0, 2000).Select(i => $"Raw file {i}\tExperiment {i % 7}\tFraction {i % 12}\tMS/MS Identified {i * 37 % 1000}")));

    private static readonly byte[] Parameters = Encoding.UTF8.GetBytes("Version\t1.6.0.16\nFixed modifications\tCarbamidomethyl (C)\n");

    /// <summary>Incompressible bytes, so the zip outgrows the 65,577-byte tail and the picked entries lie before it.</summary>
    private static readonly byte[] Evidence = RandomBytes(300_000, seed: 1609);

    private static byte[] RandomBytes(int count, int seed)
    {
        var bytes = new byte[count];
        new Random(seed).NextBytes(bytes);
        return bytes;
    }

    /// <summary>A zip built by .NET's own writer, seekable (sizes in the local headers) or streamed (data descriptors).</summary>
    private static byte[] Zip(bool streamed, params (string Name, byte[] Data, CompressionLevel Level)[] entries)
    {
        var memory = new MemoryStream();
        using (var zip = new ZipArchive(streamed ? new WriteOnlyStream(memory) : memory, ZipArchiveMode.Create, leaveOpen: true))
        {
            foreach ((string name, byte[] data, CompressionLevel level) in entries)
            {
                ZipArchiveEntry entry = zip.CreateEntry(name, level);
                if (data.Length > 0)
                {
                    using Stream stream = entry.Open();
                    stream.Write(data);
                }
            }
        }
        return memory.ToArray();
    }

    /// <summary>MaxQuant's output as a deposit zips it: a big table first, then the small design files, and a folder entry.</summary>
    private static byte[] CombinedZip(bool streamed = false) => Zip(streamed,
        ("combined/", Array.Empty<byte>(), CompressionLevel.Optimal),
        ("combined/txt/evidence.txt", Evidence, CompressionLevel.Optimal),
        ("combined/txt/summary.txt", Summary, CompressionLevel.Optimal),
        ("combined/txt/parameters.txt", Parameters, CompressionLevel.NoCompression));

    private static bool PickDesignFiles(RemoteZipEntry e) => e.Name is "summary.txt" or "parameters.txt";

    /// <summary>The offset of the central directory header of the entry named <paramref name="name"/>.</summary>
    private static int CentralHeaderOf(byte[] zip, string name)
    {
        byte[] wanted = Encoding.UTF8.GetBytes(name);
        for (int i = zip.Length - 46; i >= 0; i--)
        {
            if (BinaryPrimitives.ReadUInt32LittleEndian(zip.AsSpan(i)) == 0x02014b50
                && BinaryPrimitives.ReadUInt16LittleEndian(zip.AsSpan(i + 28)) == wanted.Length
                && zip.AsSpan(i + 46, wanted.Length).SequenceEqual(wanted))
                return i;
        }
        throw new InvalidOperationException($"No central directory header for '{name}'.");
    }

    private static int EndRecordOf(byte[] zip)
    {
        for (int i = zip.Length - 22; i >= 0; i--)
            if (BinaryPrimitives.ReadUInt32LittleEndian(zip.AsSpan(i)) == 0x06054b50)
                return i;
        throw new InvalidOperationException("No end record.");
    }

    private static void Patch16(byte[] zip, int at, int value) => BinaryPrimitives.WriteUInt16LittleEndian(zip.AsSpan(at), (ushort)value);

    private static void Patch32(byte[] zip, int at, uint value) => BinaryPrimitives.WriteUInt32LittleEndian(zip.AsSpan(at), value);

    private const long FarOffset = 5_000_000_000;
    private static readonly byte[] FarData = Encoding.ASCII.GetBytes("an entry whose local header starts beyond 4 GiB\n");

    /// <summary>
    /// A ZIP64 zip of 5 GB, sparse: one stored entry whose local header is at <see cref="FarOffset"/>, so its offset
    /// is in the ZIP64 extra field, and an end record whose fields all defer to the ZIP64 end record. Built by hand,
    /// because .NET's writer makes ZIP64 only for an entry or an offset that is really that large.
    /// </summary>
    private static RangeServer Zip64Server()
    {
        byte[] name = Encoding.ASCII.GetBytes("far/entry.txt");
        uint crc = System.IO.Hashing.Crc32.HashToUInt32(FarData);
        var bytes = new MemoryStream();
        var w = new BinaryWriter(bytes);

        w.Write(0x04034b50u); w.Write((ushort)45); w.Write((ushort)0); w.Write((ushort)0); w.Write((ushort)0); w.Write((ushort)0);
        w.Write(crc); w.Write((uint)FarData.Length); w.Write((uint)FarData.Length); w.Write((ushort)name.Length); w.Write((ushort)0);
        w.Write(name); w.Write(FarData);

        long directoryOffset = FarOffset + bytes.Length;
        long directoryStart = bytes.Length;
        w.Write(0x02014b50u); w.Write((ushort)45); w.Write((ushort)45); w.Write((ushort)0); w.Write((ushort)0); w.Write((ushort)0); w.Write((ushort)0);
        w.Write(crc); w.Write((uint)FarData.Length); w.Write((uint)FarData.Length); w.Write((ushort)name.Length); w.Write((ushort)12); w.Write((ushort)0);
        w.Write((ushort)0); w.Write((ushort)0); w.Write(0u); w.Write(0xFFFFFFFFu); w.Write(name);
        w.Write((ushort)1); w.Write((ushort)8); w.Write(FarOffset);
        long directorySize = bytes.Length - directoryStart;

        long zip64Record = FarOffset + bytes.Length;
        w.Write(0x06064b50u); w.Write(44UL); w.Write((ushort)45); w.Write((ushort)45); w.Write(0u); w.Write(0u);
        w.Write(1UL); w.Write(1UL); w.Write((ulong)directorySize); w.Write((ulong)directoryOffset);
        w.Write(0x07064b50u); w.Write(0u); w.Write((ulong)zip64Record); w.Write(1u);
        w.Write(0x06054b50u); w.Write((ushort)0xFFFF); w.Write((ushort)0xFFFF); w.Write((ushort)0xFFFF); w.Write((ushort)0xFFFF);
        w.Write(0xFFFFFFFFu); w.Write(0xFFFFFFFFu); w.Write((ushort)0);
        w.Flush();

        return new RangeServer(FarOffset + bytes.Length, (FarOffset, bytes.ToArray()));
    }

    /// <summary>A verbatim PRIDE tail (see PrideZipFixtures/PROVENANCE.md), served at its real offset in a file of its real length.</summary>
    private static RangeServer RealTailServer(string fixture, long total, string etag)
    {
        byte[] tail = File.ReadAllBytes(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "PrideZipFixtures", fixture));
        return new RangeServer(total, (total - tail.Length, tail)) { ETag = etag };
    }

    // ---- set-up -----------------------------------------------------------------

    private List<TimeSpan> _delays;
    private string _tempDir;

    /// <summary>A client over <paramref name="handler"/> whose backoffs are recorded, not slept through.</summary>
    private PrideArchiveClient ClientOver(HttpMessageHandler handler, int maxRetries = 3) =>
        new(new HttpClient(handler))
        {
            MaxRetries = maxRetries,
            RetryDelay = (delay, _) =>
            {
                _delays.Add(delay);
                return Task.CompletedTask;
            },
        };

    private static PrideArchiveFile PrideZip(string url = ZipUrl, string fileName = "combined.zip") => new()
    {
        FileName = fileName,
        FileCategory = new CvParam("PRIDE", "PRIDE:0000409", "category", "OTHER"),
        PublicFileLocations = { new CvParam("PRIDE", PrideArchiveExtensions.FtpLocationAccession, "FTP Protocol", url) },
    };

    private string Written(string relative) => Path.Combine(_tempDir, relative.Replace('/', Path.DirectorySeparatorChar));

    [SetUp]
    public void SetUp()
    {
        _delays = new List<TimeSpan>();
        _tempDir = Path.Combine(Path.GetTempPath(), "PrideZipEntriesTests", Guid.NewGuid().ToString("N"));
    }

    [TearDown]
    public void TearDown()
    {
        try
        {
            if (Directory.Exists(_tempDir))
                Directory.Delete(_tempDir, recursive: true);
        }
        catch (IOException) { }
        catch (UnauthorizedAccessException) { }
    }

    // ---- listing ------------------------------------------------------------------

    [Test]
    public async Task GetZipEntriesAsync_RealZip64TailFromPride_ListsAll41EntriesInOneRequest()
    {
        // PXD009815's mzdb_files.zip, 22.8 GB: 32 of its 41 entries start beyond 4 GiB, so their offsets exist only
        // in the ZIP64 extra field. Expected values were read from the same bytes by an independent parser.
        var server = RealTailServer("mzdb_files.zip.tail", 22_821_363_261, "\"550425e3d-5d2ba86de0300\"");
        using var client = ClientOver(server);

        IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(PrideZip(fileName: "mzdb_files.zip"));

        Assert.That(server.Requests, Is.EqualTo(new[] { ("bytes=-65577", (string)null) }));
        Assert.That(entries, Has.Count.EqualTo(41));
        Assert.That(entries[0].FullName, Is.EqualTo("mzdb_files/"));
        Assert.That(entries[0].IsDirectory, Is.True);
        Assert.That(entries[0].CanDownload, Is.False);

        RemoteZipEntry first = entries[1];
        Assert.That(first.FullName, Is.EqualTo("mzdb_files/QEKAC160601_61.raw.mzDB"));
        Assert.That(first.Name, Is.EqualTo("QEKAC160601_61.raw.mzDB"));
        Assert.That(first.CompressionMethod, Is.EqualTo(RemoteZipEntry.Deflated));
        Assert.That(first.Crc32, Is.EqualTo(0x4a9494b7u));
        Assert.That(first.CompressedLength, Is.EqualTo(591_279_325));
        Assert.That(first.Length, Is.EqualTo(1_009_557_504));
        Assert.That(first.LocalHeaderOffset, Is.EqualTo(69));

        RemoteZipEntry last = entries[40];
        Assert.That(last.FullName, Is.EqualTo("mzdb_files/QEKAC160601_46.raw.mzDB"));
        Assert.That(last.Crc32, Is.EqualTo(0x7d6cd8abu));
        Assert.That(last.CompressedLength, Is.EqualTo(578_196_859));
        Assert.That(last.LocalHeaderOffset, Is.EqualTo(22_243_161_577));
        Assert.That(entries.Count(e => e.LocalHeaderOffset > uint.MaxValue), Is.EqualTo(32));
        Assert.That(entries.Sum(e => e.CompressedLength), Is.EqualTo(22_821_354_769));
    }

    [Test]
    public async Task GetZipEntriesAsync_RealMaxQuantTailFromPride_ListsItsSixTablesInOrder()
    {
        var server = RealTailServer("MaxQuant_Results_50f-25f.zip.tail", 25_399_698, "\"1839192-5d2ba967bd080\"");
        using var client = ClientOver(server);

        IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(PrideZip(fileName: "MaxQuant_Results_50f-25f.zip"));

        Assert.That(server.Requests, Has.Count.EqualTo(1));
        Assert.That(entries.Select(e => e.FullName), Is.EqualTo(new[]
        {
            "evidence.txt", "modificationSpecificPeptides.txt", "mqpar.xml", "parameters.txt", "peptides.txt", "proteinGroups.txt",
        }));
        RemoteZipEntry mqpar = entries[2];
        Assert.That((mqpar.CompressedLength, mqpar.Length, mqpar.Crc32, mqpar.LocalHeaderOffset), Is.EqualTo((3122L, 13252L, 0xa2449417u, 20_352_315L)));
        Assert.That(entries.All(e => e.CanDownload), Is.True);
    }

    [Test]
    public async Task GetZipEntriesAsync_ZipShorterThanTheTail_IsListedFromTheWholeFileInOneRequest()
    {
        byte[] zip = Zip(false, ("a.txt", Summary, CompressionLevel.Optimal), ("b/c.txt", Parameters, CompressionLevel.NoCompression));
        Assert.That(zip.Length, Is.LessThan(65_577), "the fixture must be shorter than the tail");
        var server = new RangeServer(zip);
        using var client = ClientOver(server);

        IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(PrideZip());

        Assert.That(server.Requests, Has.Count.EqualTo(1));
        Assert.That(entries.Select(e => (e.FullName, e.CompressionMethod, e.Length)),
            Is.EqualTo(new[] { ("a.txt", RemoteZipEntry.Deflated, (long)Summary.Length), ("b/c.txt", RemoteZipEntry.Stored, (long)Parameters.Length) }));
        Assert.That(entries[1].Crc32, Is.EqualTo(System.IO.Hashing.Crc32.HashToUInt32(Parameters)));
    }

    [Test]
    public async Task GetZipEntriesAsync_DirectoryLargerThanTheTail_ReadsItWithASecondRangeUnderIfRange()
    {
        // 3,000 entries with long names make a central directory of about 280 KB, which starts before the tail.
        var names = Enumerable.Range(0, 3000).Select(i => $"combined/andromeda/a-rather-long-folder-name/run_{i:D5}.apl").ToArray();
        byte[] zip = Zip(false, names.Select(n => (n, Array.Empty<byte>(), CompressionLevel.NoCompression)).ToArray());
        int end = EndRecordOf(zip);
        long directoryOffset = BinaryPrimitives.ReadUInt32LittleEndian(zip.AsSpan(end + 16));
        long directorySize = BinaryPrimitives.ReadUInt32LittleEndian(zip.AsSpan(end + 12));
        Assert.That(directoryOffset, Is.LessThan(zip.Length - 65_577), "the fixture must put the directory before the tail");
        var server = new RangeServer(zip);
        using var client = ClientOver(server);

        IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(PrideZip());

        Assert.That(entries.Select(e => e.FullName), Is.EqualTo(names));
        Assert.That(server.Requests, Is.EqualTo(new[]
        {
            ("bytes=-65577", (string)null),
            ($"bytes={directoryOffset}-{directoryOffset + directorySize - 1}", "\"zip-v1\""),
        }));
    }

    [Test]
    public async Task GetZipEntriesAsync_Zip64WithAnEntryBeyond4GiB_ReadsItsOffsetFromTheZip64ExtraField()
    {
        var server = Zip64Server();
        using var client = ClientOver(server);

        IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(PrideZip());

        Assert.That(server.Requests, Has.Count.EqualTo(1));
        Assert.That(entries, Has.Count.EqualTo(1));
        Assert.That(entries[0].FullName, Is.EqualTo("far/entry.txt"));
        Assert.That(entries[0].LocalHeaderOffset, Is.EqualTo(FarOffset));
        Assert.That(entries[0].Length, Is.EqualTo(FarData.Length));
    }

    [Test]
    public async Task GetZipEntriesAsync_EmptyZip_ReturnsAnEmptyList()
    {
        using var client = ClientOver(new RangeServer(Zip(false)));
        Assert.That(await client.GetZipEntriesAsync(PrideZip()), Is.Empty);
    }

    [Test]
    public void GetZipEntriesAsync_NotAZip_ThrowsMzLibExceptionAfterOneRequest()
    {
        var server = new RangeServer(Evidence);
        using var client = ClientOver(server);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("'combined.zip' from ftp.pride.ebi.ac.uk is not a zip"));
        Assert.That(server.Requests, Has.Count.EqualTo(1));
        Assert.That(_delays, Is.Empty);
    }

    [Test]
    public async Task GetZipEntriesAsync_ACommentEndsTheFileWithASignatureInside_FindsTheRealEndRecord()
    {
        // A comment that itself holds the end record's signature must not be taken for it: only a record whose
        // comment runs exactly to the end of the file is.
        byte[] zip = Zip(false, ("a.txt", Parameters, CompressionLevel.Optimal));
        byte[] comment = new byte[] { 0x50, 0x4b, 0x05, 0x06 }.Concat(new byte[40]).ToArray();
        int end = EndRecordOf(zip);
        Patch16(zip, end + 20, comment.Length);
        var server = new RangeServer(zip.Concat(comment).ToArray());
        using var client = ClientOver(server);

        Assert.That((await client.GetZipEntriesAsync(PrideZip())).Single().FullName, Is.EqualTo("a.txt"));
    }

    [Test]
    public void GetZipEntriesAsync_DirectoryPointsPastTheEndRecord_ThrowsMzLibException()
    {
        byte[] zip = Zip(false, ("a.txt", Parameters, CompressionLevel.Optimal));
        Patch32(zip, EndRecordOf(zip) + 12, 1_000_000); // directory size
        using var client = ClientOver(new RangeServer(zip));

        var e = Assert.ThrowsAsync<MzLibException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("is not a readable zip: its central directory runs past its end record"));
    }

    [Test]
    public void GetZipEntriesAsync_EntryCountDisagreesWithTheDirectory_ThrowsMzLibException()
    {
        byte[] zip = Zip(false, ("a.txt", Parameters, CompressionLevel.Optimal), ("b.txt", Parameters, CompressionLevel.Optimal));
        int end = EndRecordOf(zip);
        Patch16(zip, end + 8, 3);
        Patch16(zip, end + 10, 3);
        using var client = ClientOver(new RangeServer(zip));

        var e = Assert.ThrowsAsync<MzLibException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("holds 2 entries where its end record counts 3"));
    }

    [Test]
    public void GetZipEntriesAsync_EntryOffsetPastTheData_ThrowsMzLibException()
    {
        byte[] zip = Zip(false, ("a.txt", Parameters, CompressionLevel.Optimal));
        Patch32(zip, CentralHeaderOf(zip, "a.txt") + 42, 1_000_000); // local header offset
        using var client = ClientOver(new RangeServer(zip));

        var e = Assert.ThrowsAsync<MzLibException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("'a.txt' lies outside the zip's data"));
    }

    [Test]
    public void GetZipEntriesAsync_OnePartOfASplitZip_ThrowsNotSupportedException()
    {
        byte[] zip = Zip(false, ("a.txt", Parameters, CompressionLevel.Optimal));
        Patch16(zip, EndRecordOf(zip) + 4, 1); // this disk's number
        using var client = ClientOver(new RangeServer(zip));

        var e = Assert.ThrowsAsync<NotSupportedException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("split across several files"));
    }

    [Test]
    public void GetZipEntriesAsync_ServerIgnoresRanges_ThrowsMzLibExceptionWithoutReadingTheBody()
    {
        var server = new RangeServer(CombinedZip())
        {
            Before = (_, _) => new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new UnreadableStream()) },
        };
        using var client = ClientOver(server);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("with the whole file: the server ignores byte ranges"));
        Assert.That(server.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public async Task GetZipEntriesAsync_503ThenTheTail_RetriesAndLists()
    {
        var server = new RangeServer(CombinedZip())
        {
            Before = (_, n) => n == 1 ? new HttpResponseMessage(HttpStatusCode.ServiceUnavailable) { Content = new StringContent("") } : null,
        };
        using var client = ClientOver(server);

        Assert.That(await client.GetZipEntriesAsync(PrideZip()), Has.Count.EqualTo(4));
        Assert.That(_delays, Is.EqualTo(new[] { TimeSpan.FromSeconds(5) }));
    }

    [Test]
    public void GetZipEntriesAsync_404_ThrowsHttpRequestExceptionThatSaysFailedWithStatus()
    {
        var server = new RangeServer(CombinedZip())
        {
            Before = (_, _) => new HttpResponseMessage(HttpStatusCode.NotFound) { ReasonPhrase = "Not Found", Content = new StringContent("") },
        };
        using var client = ClientOver(server);

        var e = Assert.ThrowsAsync<HttpRequestException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.StatusCode, Is.EqualTo(HttpStatusCode.NotFound));
        // pyMzLib's ClassifyError anchors on this phrase at every PRIDE status throw.
        Assert.That(e.Message, Does.Contain("failed with status 404"));
        Assert.That(ServiceFailure.IsTransient(e), Is.False);
        Assert.That(server.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public void GetZipEntriesAsync_503OnEveryAttempt_ThrowsATransientFailureOnceTheRetriesAreSpent()
    {
        var server = new RangeServer(CombinedZip())
        {
            Before = (_, _) => new HttpResponseMessage(HttpStatusCode.ServiceUnavailable) { Content = new StringContent("") },
        };
        using var client = ClientOver(server, maxRetries: 2);

        var e = Assert.ThrowsAsync<HttpRequestException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain("failed with status 503"));
        Assert.That(ServiceFailure.IsTransient(e), Is.True);
        Assert.That(server.Requests, Has.Count.EqualTo(3));
    }

    [Test]
    public void GetZipEntriesAsync_206ForOtherBytesOnEveryAttempt_ThrowsMzLibException()
    {
        byte[] zip = CombinedZip();
        var server = new RangeServer(zip);
        server.Before = (_, _) =>
        {
            var content = new ByteArrayContent(new byte[100]);
            content.Headers.ContentRange = new ContentRangeHeaderValue(0, 99, zip.Length);
            return new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = content };
        };
        using var client = ClientOver(server, maxRetries: 1);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.GetZipEntriesAsync(PrideZip()));
        Assert.That(e.Message, Does.Contain($"answer for bytes {zip.Length - 65_577}-{zip.Length - 1} of 'combined.zip' from ftp.pride.ebi.ac.uk held bytes 0-99 instead"));
        Assert.That(server.Requests, Has.Count.EqualTo(2));
    }

    [Test]
    public void GetZipEntriesAsync_ArgumentProblems_ThrowBeforeAnyRequest()
    {
        var server = new RangeServer(CombinedZip());
        using var client = ClientOver(server);

        Assert.That(Assert.ThrowsAsync<ArgumentNullException>(() => client.GetZipEntriesAsync(null)).ParamName, Is.EqualTo("zipFile"));
        Assert.That(Assert.ThrowsAsync<ArgumentException>(() => client.GetZipEntriesAsync(PrideZip(fileName: " "))).ParamName, Is.EqualTo("zipFile"));
        var asperaOnly = new PrideArchiveFile
        {
            FileName = "combined.zip",
            PublicFileLocations = { new CvParam("PRIDE", PrideArchiveExtensions.AsperaLocationAccession, "Aspera Protocol", "prd_ascp@fasp.ebi.ac.uk:pride/combined.zip") },
        };
        Assert.ThrowsAsync<NotSupportedException>(() => client.GetZipEntriesAsync(asperaOnly));
        Assert.ThrowsAsync<OperationCanceledException>(() => client.GetZipEntriesAsync(PrideZip(), new CancellationToken(canceled: true)));
        Assert.That(server.Requests, Is.Empty);
    }

    // ---- fetching -------------------------------------------------------------------

    [Test]
    public async Task DownloadZipEntriesAsync_PicksDesignFiles_WritesThemAtTheirPathsAndReadsNothingElse()
    {
        byte[] zip = CombinedZip();
        var server = new RangeServer(zip);
        using var client = ClientOver(server);
        IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(PrideZip());
        server.Requests.Clear();

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, PickDesignFiles);

        Assert.That(paths, Is.EqualTo(new[] { Written("combined/txt/summary.txt"), Written("combined/txt/parameters.txt") }));
        Assert.That(File.ReadAllBytes(paths[0]), Is.EqualTo(Summary));
        Assert.That(File.ReadAllBytes(paths[1]), Is.EqualTo(Parameters));
        Assert.That(Directory.GetFiles(_tempDir, "*", SearchOption.AllDirectories), Has.Length.EqualTo(2), "no .partial is left, and nothing else is written");

        // The tail, then each picked entry's 30-byte local header and its compressed bytes, all under the tail's ETag.
        RemoteZipEntry summary = entries.Single(e => e.Name == "summary.txt");
        RemoteZipEntry parameters = entries.Single(e => e.Name == "parameters.txt");
        long summaryData = summary.LocalHeaderOffset + 30 + "combined/txt/summary.txt".Length;
        long parametersData = parameters.LocalHeaderOffset + 30 + "combined/txt/parameters.txt".Length;
        Assert.That(server.Requests, Is.EqualTo(new[]
        {
            ("bytes=-65577", (string)null),
            ($"bytes={summary.LocalHeaderOffset}-{summary.LocalHeaderOffset + 29}", "\"zip-v1\""),
            ($"bytes={summaryData}-{summaryData + summary.CompressedLength - 1}", "\"zip-v1\""),
            ($"bytes={parameters.LocalHeaderOffset}-{parameters.LocalHeaderOffset + 29}", "\"zip-v1\""),
            ($"bytes={parametersData}-{parametersData + parameters.CompressedLength - 1}", "\"zip-v1\""),
        }));
        Assert.That(summary.CompressedLength, Is.LessThan(Summary.Length), "summary.txt must really be deflated");
        Assert.That(parameters.CompressionMethod, Is.EqualTo(RemoteZipEntry.Stored));
    }

    [Test]
    public async Task DownloadZipEntriesAsync_ZipWithDataDescriptors_TakesSizesFromTheDirectory()
    {
        // Written to a stream that cannot seek, .NET's writer leaves the local headers' sizes at zero and puts them in a
        // data descriptor after the data (flag bit 3). Only the central directory's sizes can be trusted then.
        byte[] zip = CombinedZip(streamed: true);
        Assert.That(BinaryPrimitives.ReadUInt16LittleEndian(zip.AsSpan(CentralHeaderOf(zip, "combined/txt/summary.txt") + 8)) & 8, Is.EqualTo(8),
            "the fixture must use data descriptors");
        using var client = ClientOver(new RangeServer(zip));

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, PickDesignFiles);

        Assert.That(File.ReadAllBytes(paths[0]), Is.EqualTo(Summary));
        Assert.That(File.ReadAllBytes(paths[1]), Is.EqualTo(Parameters));
    }

    [Test]
    public async Task DownloadZipEntriesAsync_EntryBeyond4GiB_IsFetchedFromItsZip64Offset()
    {
        var server = Zip64Server();
        using var client = ClientOver(server);

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, _ => true);

        Assert.That(File.ReadAllBytes(paths.Single()), Is.EqualTo(FarData));
        Assert.That(server.Requests[1].Range, Is.EqualTo($"bytes={FarOffset}-{FarOffset + 29}"));
    }

    [Test]
    public async Task DownloadZipEntriesAsync_BackslashSeparatedName_IsWrittenInItsFolder()
    {
        byte[] zip = Zip(false, (@"txt\summary.txt", Summary, CompressionLevel.Optimal));
        using var client = ClientOver(new RangeServer(zip));

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "summary.txt");

        Assert.That(paths, Is.EqualTo(new[] { Written("txt/summary.txt") }));
        Assert.That(File.ReadAllBytes(paths[0]), Is.EqualTo(Summary));
    }

    [Test]
    public async Task DownloadZipEntriesAsync_SelectorNeverSeesAFolder()
    {
        var seen = new List<string>();
        using var client = ClientOver(new RangeServer(CombinedZip()));

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e =>
        {
            seen.Add(e.FullName);
            return false;
        });

        Assert.That(seen, Is.EqualTo(new[] { "combined/txt/evidence.txt", "combined/txt/summary.txt", "combined/txt/parameters.txt" }));
        Assert.That(paths, Is.Empty);
    }

    [Test]
    public async Task DownloadZipEntriesAsync_DamagedOnTheFirstAttempt_IsFetchedAgainAndPasses()
    {
        // EBI has served damaged bytes that a fresh request got right (2026-10-08). The CRC-32 catches it here.
        byte[] zip = CombinedZip();
        var server = new RangeServer(zip) { Damage = (n, _, bytes) => { if (n == 3) bytes[^1] ^= 0xFF; } };
        using var client = ClientOver(server);

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "parameters.txt");

        Assert.That(File.ReadAllBytes(paths.Single()), Is.EqualTo(Parameters));
        Assert.That(_delays, Has.Count.EqualTo(1));
        Assert.That(server.Requests, Has.Count.EqualTo(5), "the tail, then header and data twice");
    }

    [Test]
    public void DownloadZipEntriesAsync_CrcWrongOnEveryAttempt_ThrowsMzLibExceptionAndLeavesNoFile()
    {
        var server = new RangeServer(CombinedZip()) { Damage = (n, _, bytes) => { if (n > 1 && bytes.Length > 30) bytes[^1] ^= 0xFF; } };
        using var client = ClientOver(server, maxRetries: 2);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "parameters.txt"));

        Assert.That(e.Message, Does.Contain("'combined/txt/parameters.txt' in 'combined.zip' from ftp.pride.ebi.ac.uk has CRC-32"));
        Assert.That(e.Message, Does.Contain($"where its zip records {System.IO.Hashing.Crc32.HashToUInt32(Parameters):x8}"));
        Assert.That(_delays, Has.Count.EqualTo(2));
        Assert.That(Directory.GetFiles(_tempDir, "*", SearchOption.AllDirectories), Is.Empty);
    }

    [Test]
    public void DownloadZipEntriesAsync_DeflatedBytesThatDoNotInflate_ThrowsMzLibExceptionOnceRetriesAreSpent()
    {
        byte[] zip = CombinedZip();
        var server = new RangeServer(zip) { Damage = (n, _, bytes) => { if (n > 1 && bytes.Length > 30) Array.Fill(bytes, (byte)0xFF); } };
        using var client = ClientOver(server, maxRetries: 1);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "summary.txt"));

        Assert.That(e.Message, Does.Contain("'combined/txt/summary.txt' in 'combined.zip' from ftp.pride.ebi.ac.uk did not inflate"));
        Assert.That(Directory.GetFiles(_tempDir, "*", SearchOption.AllDirectories), Is.Empty);
    }

    [Test]
    public void DownloadZipEntriesAsync_InflatesToMoreThanTheDirectorySays_StopsAndThrows()
    {
        byte[] zip = CombinedZip();
        Patch32(zip, CentralHeaderOf(zip, "combined/txt/summary.txt") + 24, 100); // the recorded length
        using var client = ClientOver(new RangeServer(zip), maxRetries: 0);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "summary.txt"));
        Assert.That(e.Message, Does.Contain("inflated to more than the 100 bytes its zip records"));
    }

    [Test]
    public void DownloadZipEntriesAsync_NoLocalHeaderWhereTheDirectoryPointsOnEveryAttempt_ThrowsMzLibException()
    {
        var server = new RangeServer(CombinedZip()) { Damage = (n, _, bytes) => { if (bytes.Length == 30) bytes[0] = 0; } };
        using var client = ClientOver(server, maxRetries: 1);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "summary.txt"));
        Assert.That(e.Message, Does.Contain("are not a zip local file header"));
    }

    [Test]
    public void DownloadZipEntriesAsync_ZipChangedBetweenReads_ThrowsMzLibExceptionWithoutReadingTheWholeFile()
    {
        // Under If-Range a changed file answers 200 with all of itself (measured on PRIDE, 2026-11-18). That body is
        // never read.
        var server = new RangeServer(CombinedZip());
        server.Before = (_, n) => n == 2 ? new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new UnreadableStream()) } : null;
        using var client = ClientOver(server);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, PickDesignFiles));
        Assert.That(e.Message, Does.Contain("with the whole file: the file has changed since its first part was read"));
        Assert.That(server.Requests, Has.Count.EqualTo(2));
        Assert.That(_delays, Is.Empty);
    }

    [Test]
    public void DownloadZipEntriesAsync_NoValidatorAndTheLengthChanged_ThrowsMzLibException()
    {
        // With no ETag and no Last-Modified there is no If-Range; a different total is how the change shows.
        byte[] zip = CombinedZip();
        var server = new RangeServer(zip) { ETag = null };
        server.Before = (request, n) =>
        {
            if (n != 2)
                return null;
            Assert.That(request.Headers.IfRange, Is.Null);
            var content = new ByteArrayContent(new byte[30]);
            RangeItemHeaderValue range = request.Headers.Range.Ranges.Single();
            content.Headers.ContentRange = new ContentRangeHeaderValue(range.From.Value, range.To.Value, zip.Length + 1);
            return new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = content };
        };
        using var client = ClientOver(server);

        var e = Assert.ThrowsAsync<MzLibException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, PickDesignFiles));
        Assert.That(e.Message, Does.Contain($"changed on the server while it was being read: it is now {zip.Length + 1} bytes, not {zip.Length}"));
    }

    [Test]
    public void DownloadZipEntriesAsync_BodyStalls_IsATransientFailureAndLeavesNoPartial()
    {
        var server = new RangeServer(CombinedZip());
        // Content-Range is checked before the body is read, so the stalled answer names exactly the bytes asked for.
        server.Before = (request, n) =>
        {
            if (n != 3)
                return null;
            RangeItemHeaderValue range = request.Headers.Range.Ranges.Single();
            var content = new StreamContent(new StallingStream());
            content.Headers.ContentRange = new ContentRangeHeaderValue(range.From.Value, range.To.Value, server.Total);
            return new HttpResponseMessage(HttpStatusCode.PartialContent) { Content = content };
        };
        using var client = new PrideArchiveClient(new HttpClient(server)) { MaxRetries = 0, BodyStallTimeout = TimeSpan.FromMilliseconds(200) };

        var e = Assert.ThrowsAsync<HttpRequestException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "summary.txt"));
        Assert.That(e.Message, Does.Contain("delivered nothing for"));
        Assert.That(ServiceFailure.IsTransient(e), Is.True);
        Assert.That(Directory.GetFiles(_tempDir, "*", SearchOption.AllDirectories), Is.Empty);
    }

    [TestCase("../evil.txt")]
    [TestCase("txt/../../evil.txt")]
    [TestCase("/rooted.txt")]
    [TestCase("C:/drive.txt")]
    [TestCase("txt//empty-segment.txt")]
    public void DownloadZipEntriesAsync_PathThatWouldLeaveTheDestination_IsRefusedBeforeAnyEntryIsRead(string name)
    {
        byte[] zip = Zip(false, ("safe.txt", Parameters, CompressionLevel.Optimal), (name, Parameters, CompressionLevel.Optimal));
        var server = new RangeServer(zip);
        using var client = ClientOver(server);

        var e = Assert.ThrowsAsync<ArgumentException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, _ => true));
        Assert.That(e.ParamName, Is.EqualTo("selector"));
        Assert.That(e.Message, Does.Contain($"'{name}' in 'combined.zip'"));
        Assert.That(server.Requests, Has.Count.EqualTo(1), "only the directory was read; not even safe.txt was fetched");
        Assert.That(Directory.Exists(_tempDir), Is.False);
    }

    [TestCase(9, false, "compressed by zip method 9")]
    [TestCase(8, true, "encrypted")]
    public async Task DownloadZipEntriesAsync_EntryThatCannotBeFetched_IsListedButRefusedBeforeAnyEntryIsRead(int method, bool encrypted, string why)
    {
        byte[] zip = CombinedZip();
        int header = CentralHeaderOf(zip, "combined/txt/summary.txt");
        Patch16(zip, header + 10, method);
        if (encrypted)
            Patch16(zip, header + 8, BinaryPrimitives.ReadUInt16LittleEndian(zip.AsSpan(header + 8)) | 1);
        var server = new RangeServer(zip);
        using var client = ClientOver(server);

        RemoteZipEntry listed = (await client.GetZipEntriesAsync(PrideZip())).Single(e => e.Name == "summary.txt");
        Assert.That((listed.CompressionMethod, listed.IsEncrypted, listed.CanDownload), Is.EqualTo((method, encrypted, false)));
        server.Requests.Clear();

        var e = Assert.ThrowsAsync<NotSupportedException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, PickDesignFiles));
        Assert.That(e.Message, Does.Contain($"'combined/txt/summary.txt' in 'combined.zip' from ftp.pride.ebi.ac.uk is {why}"));
        Assert.That(server.Requests, Has.Count.EqualTo(1));
    }

    [Test]
    public async Task DownloadZipEntriesAsync_OverwriteFalse_KeepsAnExistingFileAndDoesNotFetchIt()
    {
        string existing = Written("combined/txt/summary.txt");
        Directory.CreateDirectory(Path.GetDirectoryName(existing));
        File.WriteAllText(existing, "kept");
        var server = new RangeServer(CombinedZip());
        using var client = ClientOver(server);

        IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, PickDesignFiles, overwrite: false);

        Assert.That(paths, Is.EqualTo(new[] { existing, Written("combined/txt/parameters.txt") }));
        Assert.That(File.ReadAllText(existing), Is.EqualTo("kept"));
        Assert.That(server.Requests, Has.Count.EqualTo(3), "the tail, then parameters.txt's header and data only");
    }

    [Test]
    public async Task DownloadZipEntriesAsync_OverwriteTrue_ReplacesAnExistingFile()
    {
        string existing = Written("combined/txt/summary.txt");
        Directory.CreateDirectory(Path.GetDirectoryName(existing));
        File.WriteAllText(existing, "stale");
        using var client = ClientOver(new RangeServer(CombinedZip()));

        await client.DownloadZipEntriesAsync(PrideZip(), _tempDir, e => e.Name == "summary.txt");

        Assert.That(File.ReadAllBytes(existing), Is.EqualTo(Summary));
    }

    [Test]
    public void DownloadZipEntriesAsync_ArgumentProblems_ThrowBeforeAnyRequest()
    {
        var server = new RangeServer(CombinedZip());
        using var client = ClientOver(server);

        Assert.That(Assert.ThrowsAsync<ArgumentNullException>(() => client.DownloadZipEntriesAsync(null, _tempDir, _ => true)).ParamName, Is.EqualTo("zipFile"));
        Assert.That(Assert.ThrowsAsync<ArgumentException>(() => client.DownloadZipEntriesAsync(PrideZip(), " ", _ => true)).ParamName, Is.EqualTo("destinationDirectory"));
        Assert.That(Assert.ThrowsAsync<ArgumentNullException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, null)).ParamName, Is.EqualTo("selector"));
        Assert.ThrowsAsync<OperationCanceledException>(() => client.DownloadZipEntriesAsync(PrideZip(), _tempDir, _ => true, cancellationToken: new CancellationToken(canceled: true)));
        Assert.That(server.Requests, Is.Empty);
    }

    [Test]
    public void ZipFailures_NeverPutTheUrlInAMessage()
    {
        // A reviewer-token href carries the token in its query string. No message may carry the URL.
        const string secretUrl = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD000001/combined.zip?token=S3CR3T-REVIEWER-TOKEN";
        PrideArchiveFile file = PrideZip(secretUrl);
        byte[] zip = CombinedZip();
        var failures = new List<Exception>();

        void Capture(RangeServer server, Func<PrideArchiveClient, Task> act, int maxRetries = 1)
        {
            using var client = ClientOver(server, maxRetries);
            failures.Add(Assert.CatchAsync(() => act(client)));
        }

        Capture(new RangeServer(zip) { Before = (_, _) => new HttpResponseMessage(HttpStatusCode.Forbidden) { Content = new StringContent("") } },
            c => c.GetZipEntriesAsync(file));
        Capture(new RangeServer(Evidence), c => c.GetZipEntriesAsync(file));
        Capture(new RangeServer(zip) { Before = (_, n) => n == 2 ? new HttpResponseMessage(HttpStatusCode.OK) { Content = new StreamContent(new UnreadableStream()) } : null },
            c => c.DownloadZipEntriesAsync(file, _tempDir, PickDesignFiles));
        Capture(new RangeServer(zip) { Damage = (n, _, bytes) => { if (n > 1 && bytes.Length > 30) bytes[^1] ^= 0xFF; } },
            c => c.DownloadZipEntriesAsync(file, _tempDir, e => e.Name == "parameters.txt"));
        Capture(new RangeServer(Zip(false, ("../evil.txt", Parameters, CompressionLevel.Optimal))),
            c => c.DownloadZipEntriesAsync(file, _tempDir, _ => true));

        Assert.That(failures, Has.Count.EqualTo(5));
        foreach (Exception failure in failures)
        {
            Assert.That(failure.ToString(), Does.Not.Contain("S3CR3T"));
            Assert.That(failure.ToString(), Does.Not.Contain("https://"));
            Assert.That(failure.Message, Does.Contain("from ftp.pride.ebi.ac.uk"));
        }
    }
}

/// <summary>
/// Live canaries for listing and fetching zip entries by range against PRIDE's FTP host over HTTPS. Carries
/// [Category("ExternalService")] so CI runs it in the non-blocking external-service job, and [Category("Pride")] to
/// select it alone. Routed through <see cref="ExternalServiceTestHelper.RunAsync"/> so an EBI outage SKIPS, while a
/// contract break (not a zip, a whole-file answer to a range, wrong bytes) FAILS. Both read PXD009815 (public since
/// 2020), measured on 2026-11-18.
/// </summary>
[TestFixture]
[Category("ExternalService")]
[Category("Pride")]
[ExcludeFromCodeCoverage]
public class PrideZipEntriesLiveTests
{
    private static PrideArchiveFile Pxd009815(string fileName) => new()
    {
        FileName = fileName,
        PublicFileLocations =
        {
            new CvParam("PRIDE", PrideArchiveExtensions.FtpLocationAccession, "FTP Protocol",
                $"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2020/03/PXD009815/{fileName}"),
        },
    };

    [Test]
    public Task GetZipEntriesAsync_Live22GBZip64_ListsItsEntriesFromItsLastBytes() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            // one attempt: a PRIDE outage skips at once, not after the 5/20/60 s retry backoffs
            using var client = new PrideArchiveClient { MaxRetries = 0 };

            IReadOnlyList<RemoteZipEntry> entries = await client.GetZipEntriesAsync(Pxd009815("mzdb_files.zip"));

            Assert.That(entries, Has.Count.EqualTo(41));
            RemoteZipEntry first = entries.Single(e => e.FullName == "mzdb_files/QEKAC160601_61.raw.mzDB");
            Assert.That((first.CompressedLength, first.Length, first.Crc32), Is.EqualTo((591_279_325L, 1_009_557_504L, 0x4a9494b7u)));
            Assert.That(entries.Count(e => e.LocalHeaderOffset > uint.MaxValue), Is.EqualTo(32));
        });

    [Test]
    public Task DownloadZipEntriesAsync_LiveMaxQuantZip_FetchesMqparAndParameters() =>
        ExternalServiceTestHelper.RunAsync("PRIDE", async () =>
        {
            using var client = new PrideArchiveClient { MaxRetries = 0 };
            string dir = Path.Combine(Path.GetTempPath(), "PrideZipEntriesLive", Guid.NewGuid().ToString("N"));
            try
            {
                IReadOnlyList<string> paths = await client.DownloadZipEntriesAsync(Pxd009815("MaxQuant_Results_50f-25f.zip"), dir,
                    e => e.Name is "mqpar.xml" or "parameters.txt");

                Assert.That(paths.Select(Path.GetFileName), Is.EqualTo(new[] { "mqpar.xml", "parameters.txt" }));
                Assert.That(new FileInfo(paths[0]).Length, Is.EqualTo(13_252));
                Assert.That(new FileInfo(paths[1]).Length, Is.EqualTo(3_384));
                Assert.That(File.ReadAllText(paths[0]), Does.Contain("<MaxQuantParams"));
            }
            finally
            {
                if (Directory.Exists(dir))
                    Directory.Delete(dir, recursive: true);
            }
        });
}
