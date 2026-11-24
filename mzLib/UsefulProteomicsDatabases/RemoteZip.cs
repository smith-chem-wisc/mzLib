using System;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The zip format, as far as reading a zip on a server by byte ranges needs it: its end record (ZIP64 included),
    /// its central directory, and an entry's local header. It parses bytes and knows nothing of HTTP, so any client
    /// that can read a range of a file can list and fetch the entries of a zip with it.
    /// </summary>
    /// <remarks>
    /// The layout is PKWARE's APPNOTE.TXT (sections 4.3.7, 4.3.12, 4.3.14-4.3.16 and 4.5.3). A zip is read from its end:
    /// the end record, and before it the ZIP64 locator when there is one, sit in the last
    /// <see cref="TailLength"/> bytes, and they give the central directory's offset and size. In the zips measured on
    /// PRIDE on 2026-11-18 (PXD009815's 22.8 GB ZIP64 <c>mzdb_files.zip</c>, 41 entries, and its 25 MB MaxQuant
    /// zip, 6 entries) the whole directory sat inside those bytes too, so one read listed them.
    /// <para>
    /// A zip whose end record cannot be found, or whose directory points outside the zip, is an
    /// <see cref="MzLibException"/>: the server answered, and the bytes are not a readable zip. A zip split across
    /// several files is a <see cref="NotSupportedException"/>.
    /// </para>
    /// </remarks>
    internal static class RemoteZip
    {
        /// <summary>The end-of-central-directory record's length, without its comment.</summary>
        internal const int EndRecordLength = 22;

        /// <summary>The ZIP64 end-of-central-directory locator's length. It sits just before the end record.</summary>
        internal const int Zip64LocatorLength = 20;

        /// <summary>The ZIP64 end-of-central-directory record's length, without its extensible data.</summary>
        internal const int Zip64EndRecordLength = 56;

        /// <summary>A local file header's length, without its name and extra field.</summary>
        internal const int LocalHeaderLength = 30;

        /// <summary>
        /// How many bytes from a zip's end are read first: the end record with the longest comment it can carry
        /// (65,535 bytes), and the ZIP64 locator before it. An end record not found in them means the file is not a zip.
        /// </summary>
        internal const int TailLength = Zip64LocatorLength + EndRecordLength + ushort.MaxValue;

        private const int CentralHeaderLength = 46;
        private const uint EndRecordSignature = 0x06054b50;
        private const uint Zip64LocatorSignature = 0x07064b50;
        private const uint Zip64EndRecordSignature = 0x06064b50;
        private const uint CentralHeaderSignature = 0x02014b50;
        private const uint LocalHeaderSignature = 0x04034b50;
        private const ushort Zip64ExtraId = 0x0001;
        private const uint Unknown32 = 0xFFFFFFFF;
        private const ushort Unknown16 = 0xFFFF;

        /// <summary>Where a zip's central directory is, as its end records say.</summary>
        /// <param name="DirectoryOffset">The directory's first byte.</param>
        /// <param name="DirectorySize">The directory's length in bytes.</param>
        /// <param name="EntryCount">How many entries the end record counts.</param>
        /// <param name="EndOffset">Where the record that bounds the directory starts (the ZIP64 record, else the end record).</param>
        /// <param name="Zip64RecordOffset">Where the ZIP64 end record is, while it has still to be read; else null.</param>
        internal readonly record struct ZipEnd(long DirectoryOffset, long DirectorySize, long EntryCount, long EndOffset,
            long? Zip64RecordOffset);

        /// <summary>
        /// Finds the end record in the last bytes of a zip. When the zip is ZIP64 the result names the ZIP64 end record
        /// still to be read (<see cref="ZipEnd.Zip64RecordOffset"/>), for <see cref="ReadZip64EndRecord"/>.
        /// </summary>
        /// <param name="tail">The zip's last bytes, ending at its last byte.</param>
        /// <param name="tailStart">The offset of <paramref name="tail"/>'s first byte in the zip.</param>
        /// <param name="described">How messages name the zip.</param>
        internal static ZipEnd ReadEndRecord(ReadOnlySpan<byte> tail, long tailStart, string described)
        {
            // From the end backwards, the first signature whose comment runs exactly to the end of the file. Requiring
            // that keeps a signature that happens to lie inside the comment, or inside compressed data, from being taken.
            int first = Math.Max(0, tail.Length - EndRecordLength - ushort.MaxValue);
            for (int i = tail.Length - EndRecordLength; i >= first; i--)
            {
                if (U32(tail, i) != EndRecordSignature || i + EndRecordLength + U16(tail, i + 20) != tail.Length)
                    continue;

                long at = tailStart + i;
                if (i >= Zip64LocatorLength && U32(tail, i - Zip64LocatorLength) == Zip64LocatorSignature)
                {
                    int locator = i - Zip64LocatorLength;
                    if (U32(tail, locator + 4) != 0 || U32(tail, locator + 16) > 1)
                        throw Split(described);
                    ulong record = U64(tail, locator + 8);
                    if (record > (ulong)(at - Zip64LocatorLength - Zip64EndRecordLength))
                        throw Damaged(described, "its ZIP64 locator points past itself");
                    return new ZipEnd(0, 0, 0, (long)record, (long)record);
                }

                long count = U16(tail, i + 10);
                long size = U32(tail, i + 12);
                long offset = U32(tail, i + 16);
                if (U16(tail, i + 4) != 0 || U16(tail, i + 6) != 0 || U16(tail, i + 8) != count)
                    throw Split(described);
                if (offset + size > at)
                    throw Damaged(described, "its central directory runs past its end record");
                return new ZipEnd(offset, size, count, at, null);
            }

            throw new MzLibException($"{described} is not a zip: its last {tail.Length} bytes hold no end-of-central-directory record.");
        }

        /// <summary>Completes a ZIP64 <paramref name="end"/> from the ZIP64 end record its locator points to.</summary>
        /// <param name="record">The bytes at <see cref="ZipEnd.Zip64RecordOffset"/>, at least <see cref="Zip64EndRecordLength"/> of them.</param>
        /// <param name="end">What <see cref="ReadEndRecord"/> found.</param>
        /// <param name="described">How messages name the zip.</param>
        internal static ZipEnd ReadZip64EndRecord(ReadOnlySpan<byte> record, ZipEnd end, string described)
        {
            if (record.Length < Zip64EndRecordLength || U32(record, 0) != Zip64EndRecordSignature)
                throw Damaged(described, "there is no ZIP64 end record where its locator points");
            if (U32(record, 16) != 0 || U32(record, 20) != 0 || U64(record, 24) != U64(record, 32))
                throw Split(described);

            ulong count = U64(record, 32);
            ulong size = U64(record, 40);
            ulong offset = U64(record, 48);
            ulong at = (ulong)end.EndOffset;
            if (count > long.MaxValue || size > at || offset > at - size)
                throw Damaged(described, "its central directory runs past its ZIP64 end record");
            return new ZipEnd((long)offset, (long)size, (long)count, end.EndOffset, null);
        }

        /// <summary>Reads every entry of a central directory, in the order the directory lists them.</summary>
        /// <param name="directory">The whole directory: <see cref="ZipEnd.DirectorySize"/> bytes from <see cref="ZipEnd.DirectoryOffset"/>.</param>
        /// <param name="end">Where the directory is, and how many entries the end record counts.</param>
        /// <param name="described">How messages name the zip.</param>
        internal static List<RemoteZipEntry> ReadCentralDirectory(ReadOnlySpan<byte> directory, ZipEnd end, string described)
        {
            var entries = new List<RemoteZipEntry>((int)Math.Min(end.EntryCount, 65536));
            int p = 0;
            while (p < directory.Length)
            {
                if (directory.Length - p < CentralHeaderLength || U32(directory, p) != CentralHeaderSignature)
                    throw Damaged(described, $"its central directory breaks off after {entries.Count} entries");

                int flags = U16(directory, p + 8);
                int method = U16(directory, p + 10);
                uint crc = U32(directory, p + 16);
                long compressedLength = U32(directory, p + 20);
                long length = U32(directory, p + 24);
                int nameLength = U16(directory, p + 28);
                int extraLength = U16(directory, p + 30);
                int commentLength = U16(directory, p + 32);
                long disk = U16(directory, p + 34);
                long localHeaderOffset = U32(directory, p + 42);
                if (directory.Length - p - CentralHeaderLength < nameLength + extraLength + commentLength)
                    throw Damaged(described, $"its central directory breaks off after {entries.Count} entries");

                string name = Encoding.UTF8.GetString(directory.Slice(p + CentralHeaderLength, nameLength));
                ReadOnlySpan<byte> extra = directory.Slice(p + CentralHeaderLength + nameLength, extraLength);

                // A field too large for 32 (or 16) bits holds all ones, and its value is in the ZIP64 extra field, in
                // this order, present only for the fields that overflowed (APPNOTE 4.5.3).
                if (length == Unknown32 || compressedLength == Unknown32 || localHeaderOffset == Unknown32 || disk == Unknown16)
                {
                    ReadOnlySpan<byte> zip64 = FindExtra(extra, Zip64ExtraId);
                    int q = 0;
                    if (length == Unknown32 && TryTake64(zip64, ref q, out long value64))
                        length = value64;
                    if (compressedLength == Unknown32 && TryTake64(zip64, ref q, out value64))
                        compressedLength = value64;
                    if (localHeaderOffset == Unknown32 && TryTake64(zip64, ref q, out value64))
                        localHeaderOffset = value64;
                    if (disk == Unknown16 && zip64.Length - q >= 4)
                        disk = U32(zip64, q);
                    if (length == Unknown32 || compressedLength == Unknown32 || localHeaderOffset == Unknown32 || disk == Unknown16)
                        throw Damaged(described, $"'{name}' has no ZIP64 field for a size or offset too large for its header");
                }

                if (disk != 0)
                    throw Split(described);
                if (length < 0 || compressedLength < 0 || localHeaderOffset < 0
                    || localHeaderOffset > end.DirectoryOffset - LocalHeaderLength
                    || compressedLength > end.DirectoryOffset - LocalHeaderLength - localHeaderOffset)
                    throw Damaged(described, $"'{name}' lies outside the zip's data");

                entries.Add(new RemoteZipEntry(name, length, compressedLength, crc, method, (flags & 1) != 0, localHeaderOffset));
                p += CentralHeaderLength + nameLength + extraLength + commentLength;
            }

            // The end record's 16-bit count wraps in a zip of more than 65,535 entries that its writer did not make ZIP64.
            // The directory's size is what bounds it then, so the count is checked modulo 65,536.
            if (entries.Count != end.EntryCount && entries.Count % 65536 != end.EntryCount)
                throw Damaged(described, $"its central directory holds {entries.Count} entries where its end record counts {end.EntryCount}");
            return entries;
        }

        /// <summary>
        /// Where an entry's data starts, from its local header (which carries its own name and extra field, of their own
        /// lengths). Bytes that are not a local header are a <see cref="DownloadEngine"/> content mismatch, retried as
        /// damaged bytes are: the directory said a header is there.
        /// </summary>
        /// <param name="header">The <see cref="LocalHeaderLength"/> bytes at the entry's <see cref="RemoteZipEntry.LocalHeaderOffset"/>.</param>
        /// <param name="entry">The entry.</param>
        /// <param name="directoryOffset">Where the central directory starts: no entry's data runs past it.</param>
        /// <param name="described">How messages name the entry.</param>
        internal static long ReadLocalHeader(ReadOnlySpan<byte> header, RemoteZipEntry entry, long directoryOffset, string described)
        {
            if (header.Length < LocalHeaderLength || U32(header, 0) != LocalHeaderSignature)
                throw new ContentMismatchException($"The bytes read for {described} are not a zip local file header.");

            long dataStart = entry.LocalHeaderOffset + LocalHeaderLength + U16(header, 26) + U16(header, 28);
            if (dataStart > directoryOffset - entry.CompressedLength)
                throw Damaged(described, "its data runs into the central directory");
            return dataStart;
        }

        /// <summary>The data of the extra-field block with header <paramref name="id"/>, or empty when there is none.</summary>
        private static ReadOnlySpan<byte> FindExtra(ReadOnlySpan<byte> extra, ushort id)
        {
            int e = 0;
            while (extra.Length - e >= 4)
            {
                int size = U16(extra, e + 2);
                if (extra.Length - e - 4 < size)
                    break;
                if (U16(extra, e) == id)
                    return extra.Slice(e + 4, size);
                e += 4 + size;
            }
            return ReadOnlySpan<byte>.Empty;
        }

        private static bool TryTake64(ReadOnlySpan<byte> data, ref int position, out long value)
        {
            value = 0;
            if (data.Length - position < 8 || U64(data, position) > long.MaxValue)
                return false;
            value = (long)U64(data, position);
            position += 8;
            return true;
        }

        private static MzLibException Damaged(string described, string reason) =>
            new($"{described} is not a readable zip: {reason}.");

        private static NotSupportedException Split(string described) =>
            new($"{described} is one part of a zip split across several files, which cannot be read on its own.");

        private static ushort U16(ReadOnlySpan<byte> bytes, int at) => BinaryPrimitives.ReadUInt16LittleEndian(bytes.Slice(at));

        private static uint U32(ReadOnlySpan<byte> bytes, int at) => BinaryPrimitives.ReadUInt32LittleEndian(bytes.Slice(at));

        private static ulong U64(ReadOnlySpan<byte> bytes, int at) => BinaryPrimitives.ReadUInt64LittleEndian(bytes.Slice(at));
    }
}
