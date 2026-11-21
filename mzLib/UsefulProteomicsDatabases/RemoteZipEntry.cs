using System;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One entry of a zip file on a server, as its central directory records it: read by byte ranges, without
    /// downloading the zip. Returned by <see cref="PrideArchiveClient.GetZipEntriesAsync"/>, and handed to the selector
    /// of <see cref="PrideArchiveClient.DownloadZipEntriesAsync"/>.
    /// </summary>
    /// <remarks>
    /// The names follow <see cref="System.IO.Compression.ZipArchiveEntry"/>'s: <see cref="FullName"/>, <see cref="Name"/>,
    /// <see cref="Length"/>, <see cref="CompressedLength"/> and <see cref="Crc32"/> mean what they mean there.
    /// </remarks>
    public sealed class RemoteZipEntry
    {
        /// <summary>The compression method of an entry stored as it is (zip method 0).</summary>
        public const int Stored = 0;

        /// <summary>The compression method of a deflated entry (zip method 8), the one nearly every zip uses.</summary>
        public const int Deflated = 8;

        internal RemoteZipEntry(string fullName, long length, long compressedLength, uint crc32, int compressionMethod,
            bool isEncrypted, long localHeaderOffset)
        {
            FullName = fullName ?? throw new ArgumentNullException(nameof(fullName));
            Length = length;
            CompressedLength = compressedLength;
            Crc32 = crc32;
            CompressionMethod = compressionMethod;
            IsEncrypted = isEncrypted;
            LocalHeaderOffset = localHeaderOffset;
        }

        /// <summary>
        /// The entry's path inside the zip, as the zip records it, e.g. <c>combined/txt/summary.txt</c>. Folders are
        /// separated by <c>/</c> (some Windows tools write <c>\</c>). A folder entry ends with its separator.
        /// </summary>
        /// <remarks>
        /// Decoded as UTF-8, as <see cref="System.IO.Compression.ZipArchive"/> does on .NET when the zip does not set
        /// its UTF-8 flag. An ASCII name reads the same either way.
        /// </remarks>
        public string FullName { get; }

        /// <summary>The last segment of <see cref="FullName"/>: the file's own name, e.g. <c>summary.txt</c>. Empty for a folder.</summary>
        public string Name => FullName.Substring(FullName.LastIndexOfAny(new[] { '/', '\\' }) + 1);

        /// <summary>The entry's size once inflated, in bytes.</summary>
        public long Length { get; }

        /// <summary>The entry's size inside the zip, in bytes: what fetching it costs.</summary>
        public long CompressedLength { get; }

        /// <summary>The CRC-32 the zip records for the inflated entry. A fetched entry is checked against it.</summary>
        public uint Crc32 { get; }

        /// <summary>
        /// The zip compression method: <see cref="Stored"/> or <see cref="Deflated"/>, the two that can be fetched, or
        /// another (9 Deflate64, 12 BZIP2, 14 LZMA, 93 Zstandard ...), which is listed but cannot be.
        /// </summary>
        public int CompressionMethod { get; }

        /// <summary>Whether the entry is encrypted. An encrypted entry is listed but cannot be fetched.</summary>
        public bool IsEncrypted { get; }

        /// <summary>Whether the entry is a folder (its name ends with a separator). A folder has nothing to fetch.</summary>
        public bool IsDirectory => FullName.EndsWith('/') || FullName.EndsWith('\\');

        /// <summary>Whether <see cref="PrideArchiveClient.DownloadZipEntriesAsync"/> can fetch the entry: a file, not encrypted, stored or deflated.</summary>
        public bool CanDownload => !IsDirectory && !IsEncrypted && CompressionMethod is Stored or Deflated;

        /// <summary>Where the entry's local header starts in the zip: the first byte a fetch of it reads.</summary>
        internal long LocalHeaderOffset { get; }

        /// <inheritdoc />
        public override string ToString() => FullName;
    }
}
