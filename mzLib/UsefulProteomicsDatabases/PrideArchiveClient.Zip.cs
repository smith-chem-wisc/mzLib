using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.IO.Hashing;
using System.Linq;
using System.Net.Http;
using System.Net.Http.Headers;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    public sealed partial class PrideArchiveClient
    {
        /// <summary>
        /// Lists the entries of a zip file in a PRIDE project without downloading it: only the zip's central
        /// directory is read, by HTTP range requests on its HTTPS location.
        /// </summary>
        /// <remarks>
        /// One request reads the last 65,577 bytes of the zip (the end record with the longest comment it can carry, and
        /// the ZIP64 locator before it), and the length of the whole file comes with them. A second request is made only
        /// when the ZIP64 end record or the central directory lies before those bytes. Measured on 2026-11-18, PXD009815's
        /// 22.8 GB ZIP64 <c>mzdb_files.zip</c> (41 entries) and its 25 MB MaxQuant zip (6 entries) were both listed by
        /// the one request. Later requests carry the first one's validator (<c>If-Range</c>), so a zip PRIDE rewrites
        /// between them is reported, never read as a mix of two versions.
        /// <para>
        /// Retries, timeouts and the stall deadline are those of <see cref="DownloadFileAsync(PrideArchiveFile, string, bool, CancellationToken)"/>.
        /// </para>
        /// </remarks>
        /// <param name="zipFile">The zip, e.g. from <see cref="GetProjectFilesAsync"/>. Must not be null and must have a file name.</param>
        /// <param name="cancellationToken">Cancels the listing.</param>
        /// <returns>Every entry, folders included, in the order the central directory lists them. Never null; empty for an empty zip.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="zipFile"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="zipFile"/> has no file name.</exception>
        /// <exception cref="NotSupportedException">The file has no HTTPS-reachable location (e.g. Aspera-only), or it is one part of a zip split across several files.</exception>
        /// <exception cref="HttpRequestException">A request returned a non-success status code (carried in <see cref="HttpRequestException.StatusCode"/>; the message says "failed with status NNN"), PRIDE did not answer within the client's timeout, or a body dropped or stalled. A transient one is thrown only after its retries are spent. The message names the file and host, never the URL.</exception>
        /// <exception cref="MzLibException">The file is not a readable zip; PRIDE answered a range request with the whole file (the zip changed while it was read, or the server ignores ranges); or an answer held other bytes than those asked for, on every attempt.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled.</exception>
        public async Task<IReadOnlyList<RemoteZipEntry>> GetZipEntriesAsync(PrideArchiveFile zipFile, CancellationToken cancellationToken = default)
        {
            CheckZipFile(zipFile);
            cancellationToken.ThrowIfCancellationRequested();

            string url = zipFile.GetHttpsDownloadUrl(); // throws NotSupportedException if unreachable over HTTPS
            string described = DownloadEngine.DescribeDownload(url, zipFile.FileName);
            try
            {
                return (await ReadZipDirectoryAsync(url, described, cancellationToken).ConfigureAwait(false)).Entries;
            }
            catch (ContentMismatchException e)
            {
                throw new MzLibException(e.Message, e);
            }
        }

        /// <summary>
        /// Fetches the entries of a zip file in a PRIDE project that <paramref name="selector"/> picks, by HTTP range
        /// requests, without downloading the rest of the zip. Each is inflated, checked against the CRC-32 and the size
        /// the zip records for it, and written below <paramref name="destinationDirectory"/> at its path inside the zip.
        /// </summary>
        /// <remarks>
        /// The central directory is read as <see cref="GetZipEntriesAsync"/> reads it. Then each entry costs two requests:
        /// its 30-byte local header (which gives where its data starts), then its compressed bytes, streamed through the
        /// inflater into a ".partial" file that is moved into place only once the CRC-32 and size match. For MaxQuant's
        /// <c>mqpar.xml</c> and <c>parameters.txt</c> in PXD009815's 25 MB zip, that is 4 KB of the zip.
        /// <para>
        /// Bytes that fail the check (EBI has served damaged bytes that a fresh request got right, 2026-10-08) are fetched
        /// again from the entry's first byte, within the same retry budget as a transient failure. Only when every attempt
        /// fails does it throw <see cref="MzLibException"/>. An entry is never resumed: a failed one's ".partial" is deleted.
        /// </para>
        /// <para>
        /// Everything that can be refused is refused before any entry is fetched: an entry that cannot be fetched (see
        /// <see cref="RemoteZipEntry.CanDownload"/>), and an entry whose path could not be written below the destination.
        /// </para>
        /// </remarks>
        /// <param name="zipFile">The zip, e.g. from <see cref="GetProjectFilesAsync"/>. Must not be null and must have a file name.</param>
        /// <param name="destinationDirectory">The directory to write into; it and the entries' folders are created if they do not exist.</param>
        /// <param name="selector">
        /// Picks the entries to fetch. It sees every file entry, never a folder. For MaxQuant's design files:
        /// <c>e =&gt; e.Name is "summary.txt" or "mqpar.xml" or "parameters.txt"</c>.
        /// </param>
        /// <param name="overwrite">
        /// When true (the default) an existing file is replaced. When false, an entry whose file already exists is left
        /// untouched and not fetched; the central directory is still read.
        /// </param>
        /// <param name="cancellationToken">Cancels the listing and the fetches.</param>
        /// <returns>The full path of each picked entry's file, in the order the central directory lists them. Never null; empty when nothing was picked.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="zipFile"/> or <paramref name="selector"/> is null.</exception>
        /// <exception cref="ArgumentException"><paramref name="zipFile"/> has no file name, the destination directory is blank, or a picked entry's path has a segment that cannot be a local file name (a "..", a root, a drive, an empty segment).</exception>
        /// <exception cref="NotSupportedException">As for <see cref="GetZipEntriesAsync"/>, or a picked entry is encrypted or compressed by a method other than stored or deflated.</exception>
        /// <exception cref="HttpRequestException">As for <see cref="GetZipEntriesAsync"/>.</exception>
        /// <exception cref="MzLibException">As for <see cref="GetZipEntriesAsync"/>, or an entry's bytes did not inflate to its recorded size and CRC-32 on every attempt. The message names the entry, the zip and the host, never the URL.</exception>
        /// <exception cref="IOException">Writing a file failed (the caller's disk). Not retried.</exception>
        /// <exception cref="OperationCanceledException"><paramref name="cancellationToken"/> was cancelled.</exception>
        public async Task<IReadOnlyList<string>> DownloadZipEntriesAsync(PrideArchiveFile zipFile, string destinationDirectory,
            Func<RemoteZipEntry, bool> selector, bool overwrite = true, CancellationToken cancellationToken = default)
        {
            CheckZipFile(zipFile);
            if (string.IsNullOrWhiteSpace(destinationDirectory))
                throw new ArgumentException("A destination directory is required.", nameof(destinationDirectory));
            if (selector == null)
                throw new ArgumentNullException(nameof(selector));
            cancellationToken.ThrowIfCancellationRequested();

            string url = zipFile.GetHttpsDownloadUrl(); // throws NotSupportedException if unreachable over HTTPS
            string described = DownloadEngine.DescribeDownload(url, zipFile.FileName);
            try
            {
                ZipDirectory directory = await ReadZipDirectoryAsync(url, described, cancellationToken).ConfigureAwait(false);

                var picked = new List<(RemoteZipEntry Entry, string Path)>();
                foreach (RemoteZipEntry entry in directory.Entries.Where(e => !e.IsDirectory && selector(e)))
                {
                    if (!entry.CanDownload)
                        throw new NotSupportedException(
                            $"'{entry.FullName}' in {described} is {(entry.IsEncrypted ? "encrypted" : $"compressed by zip method {entry.CompressionMethod}")}; only stored and deflated entries can be fetched.");
                    picked.Add((entry, LocalPathOfZipEntry(entry, destinationDirectory, described)));
                }

                var paths = new List<string>(picked.Count);
                foreach ((RemoteZipEntry entry, string path) in picked)
                {
                    cancellationToken.ThrowIfCancellationRequested();
                    if (overwrite || !File.Exists(path))
                    {
                        Directory.CreateDirectory(Path.GetDirectoryName(path));
                        await DownloadZipEntryAsync(url, $"'{entry.FullName}' in {described}", directory, entry, path, cancellationToken).ConfigureAwait(false);
                    }
                    paths.Add(path);
                }
                return paths;
            }
            catch (ContentMismatchException e)
            {
                // Every attempt came back wrong, so this is no longer an outage worth waiting out.
                throw new MzLibException(e.Message, e);
            }
        }

        private static void CheckZipFile(PrideArchiveFile zipFile)
        {
            if (zipFile == null)
                throw new ArgumentNullException(nameof(zipFile));
            if (string.IsNullOrWhiteSpace(zipFile.FileName))
                throw new ArgumentException("The PRIDE file has no file name.", nameof(zipFile));
        }

        /// <summary>A zip's central directory, and what later reads of the same zip are checked against.</summary>
        /// <param name="Entries">The entries.</param>
        /// <param name="DirectoryOffset">Where the directory starts: no entry's data runs past it.</param>
        /// <param name="Total">The zip's length.</param>
        /// <param name="Validator">The first read's validator, sent as <c>If-Range</c> with every later read; null when PRIDE sent none.</param>
        private sealed record ZipDirectory(IReadOnlyList<RemoteZipEntry> Entries, long DirectoryOffset, long Total,
            RangeConditionHeaderValue Validator);

        private async Task<ZipDirectory> ReadZipDirectoryAsync(string url, string described, CancellationToken cancellationToken)
        {
            DownloadEngine.RangeBytes tail = await ReadZipRangeAsync(url, described, null, RemoteZip.TailLength, null, null, cancellationToken).ConfigureAwait(false);
            RangeConditionHeaderValue validator = tail.Validator;

            RemoteZip.ZipEnd end = RemoteZip.ReadEndRecord(tail.Bytes, tail.From, described);
            if (end.Zip64RecordOffset is long zip64Record)
            {
                byte[] record = await SliceOrReadAsync(tail, zip64Record, RemoteZip.Zip64EndRecordLength, url, described, validator, cancellationToken).ConfigureAwait(false);
                end = RemoteZip.ReadZip64EndRecord(record, end, described);
            }

            if (end.DirectorySize > int.MaxValue)
                throw new NotSupportedException($"The central directory of {described} ({end.DirectorySize} bytes) is too large to read into memory.");

            byte[] directory = await SliceOrReadAsync(tail, end.DirectoryOffset, (int)end.DirectorySize, url, described, validator, cancellationToken).ConfigureAwait(false);
            return new ZipDirectory(RemoteZip.ReadCentralDirectory(directory, end, described), end.DirectoryOffset, tail.Total, validator);
        }

        /// <summary>Bytes <paramref name="offset"/> to <paramref name="offset"/> + <paramref name="length"/> - 1 of the zip: from the tail already read when they lie in it, else read.</summary>
        private async Task<byte[]> SliceOrReadAsync(DownloadEngine.RangeBytes tail, long offset, int length, string url, string described,
            RangeConditionHeaderValue validator, CancellationToken cancellationToken)
        {
            if (length == 0)
                return Array.Empty<byte>();
            if (offset >= tail.From && offset - tail.From <= tail.Bytes.Length - length)
                return tail.Bytes.AsSpan((int)(offset - tail.From), length).ToArray();

            DownloadEngine.RangeBytes read = await ReadZipRangeAsync(url, described, offset, offset + length - 1, validator, tail.Total, cancellationToken).ConfigureAwait(false);
            return read.Bytes;
        }

        private Task<DownloadEngine.RangeBytes> ReadZipRangeAsync(string url, string described, long? from, long to,
            RangeConditionHeaderValue validator, long? total, CancellationToken cancellationToken) =>
            _engine.WithRetryAsync(
                () => _engine.ReadRangeOnceAsync(url, described, from, to, validator, total, cancellationToken),
                _engine.HostOf(url), cancellationToken);

        /// <summary>
        /// Fetches one entry into <paramref name="path"/>. One attempt reads the entry's local header and then its data,
        /// and a failed attempt starts again from the header.
        /// </summary>
        private async Task DownloadZipEntryAsync(string url, string described, ZipDirectory directory, RemoteZipEntry entry,
            string path, CancellationToken cancellationToken)
        {
            string partialPath = path + ".partial";
            var check = new ZipEntryCheck(entry);
            try
            {
                await _engine.WithRetryAsync(async () =>
                {
                    check.Reset();
                    DownloadEngine.RangeBytes header = await _engine.ReadRangeOnceAsync(url, described, entry.LocalHeaderOffset,
                        entry.LocalHeaderOffset + RemoteZip.LocalHeaderLength - 1, directory.Validator, directory.Total, cancellationToken).ConfigureAwait(false);
                    long dataStart = RemoteZip.ReadLocalHeader(header.Bytes, entry, directory.DirectoryOffset, described);

                    using (var file = new FileStream(partialPath, FileMode.Create, FileAccess.Write, FileShare.None))
                    {
                        if (entry.CompressedLength > 0)
                        {
                            using HttpResponseMessage response = await _engine.GetRangeAsync(url, described, dataStart,
                                dataStart + entry.CompressedLength - 1, directory.Validator, directory.Total, cancellationToken).ConfigureAwait(false);
                            using Stream body = await response.Content.ReadAsStreamAsync(cancellationToken).ConfigureAwait(false);
                            using Stream data = entry.CompressionMethod == RemoteZipEntry.Deflated
                                ? new DeflateStream(body, CompressionMode.Decompress)
                                : body;
                            try
                            {
                                await _engine.CopyUntilStalledAsync(data, file, described, check, cancellationToken).ConfigureAwait(false);
                            }
                            catch (InvalidDataException)
                            {
                                // Bytes the inflater refuses were damaged on the way, or are damaged in the zip; a fresh read tells which.
                                throw new ContentMismatchException($"The PRIDE download of {described} did not inflate.");
                            }
                        }
                    }

                    if (check.Mismatch is string mismatch)
                        throw new ContentMismatchException($"The PRIDE download of {described} {mismatch}.");
                    return true;
                }, _engine.HostOf(url), cancellationToken).ConfigureAwait(false);

                File.Move(partialPath, path, overwrite: true);
            }
            catch
            {
                // An entry is never resumed (it is small next to its zip, and a resumed inflate would need the inflater's
                // state), so a failed one's bytes have no use.
                DownloadEngine.DeleteQuietly(partialPath);
                throw;
            }
        }

        /// <summary>
        /// Where <paramref name="entry"/> is written below <paramref name="destinationDirectory"/>: at its path inside the
        /// zip, each segment of which must be a plain local file name, so a name in the zip cannot write outside the
        /// destination. Both <c>/</c> and <c>\</c> separate folders.
        /// </summary>
        private static string LocalPathOfZipEntry(RemoteZipEntry entry, string destinationDirectory, string described)
        {
            char[] invalid = Path.GetInvalidFileNameChars();
            string[] segments = entry.FullName.Split('/', '\\');
            foreach (string segment in segments)
            {
                if (segment.Length == 0 || segment is "." or ".." || segment.IndexOfAny(invalid) >= 0
                    || segment.Contains(':') || Path.IsPathRooted(segment))
                    throw new ArgumentException(
                        $"The path of '{entry.FullName}' in {described} has a segment that cannot be a local file name; refusing to write outside the destination directory.",
                        "selector");
            }

            string root = Path.GetFullPath(destinationDirectory);
            string path = Path.GetFullPath(Path.Combine(new[] { root }.Concat(segments).ToArray()));
            if (!path.StartsWith(Path.TrimEndingDirectorySeparator(root) + Path.DirectorySeparatorChar, StringComparison.Ordinal))
                throw new ArgumentException(
                    $"The path of '{entry.FullName}' in {described} resolves outside the destination directory.", "selector");
            return path;
        }

        /// <summary>Checks an entry's inflated bytes, as they arrive, against the size and CRC-32 its zip records.</summary>
        private sealed class ZipEntryCheck(RemoteZipEntry entry) : IDownloadInspector
        {
            private readonly Crc32 _crc = new();
            private long _length;

            public void Reset()
            {
                _crc.Reset();
                _length = 0;
            }

            // An entry is never resumed, so its bytes always arrive from the first.
            public void ContinueAt(long offset) => throw new NotSupportedException();

            public void Append(ReadOnlySpan<byte> data)
            {
                _crc.Append(data);
                _length += data.Length;
            }

            public string StopReason =>
                _length > entry.Length ? $"inflated to more than the {entry.Length} bytes its zip records" : null;

            /// <summary>Null when the bytes seen are the entry; else the end of a sentence saying how they differ.</summary>
            public string Mismatch =>
                _length != entry.Length ? $"inflated to {_length} bytes where its zip records {entry.Length}"
                : _crc.GetCurrentHashAsUInt32() != entry.Crc32 ? $"has CRC-32 {_crc.GetCurrentHashAsUInt32():x8} where its zip records {entry.Crc32:x8}"
                : null;
        }
    }
}
