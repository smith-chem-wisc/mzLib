using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Net.Security;
using System.Net.Sockets;
using System.Security.Authentication;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using FluentFTP;
using FluentFTP.Exceptions;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// The FTP operations <see cref="MassiveClient"/> needs, behind one seam so offline tests can stand in for the
    /// server. One connection, used by one operation at a time.
    /// </summary>
    /// <remarks>
    /// Every implementation fails in three ways only, which the client then maps onto mzLib's error contract:
    /// <list type="bullet">
    /// <item><see cref="MassiveFtpReplyException"/>: the server answered with a negative reply, carrying its code;</item>
    /// <item><see cref="AuthenticationException"/>: the TLS handshake or the certificate check failed;</item>
    /// <item><see cref="MassiveFtpTransportException"/>: anything else on the wire (a timeout, a dropped or refused
    /// connection, an unreadable reply), on a command or while reading a file.</item>
    /// </list>
    /// A caller's cancellation stays an <see cref="OperationCanceledException"/>.
    /// </remarks>
    internal interface IMassiveFtpConnection : IDisposable
    {
        /// <summary>The raw lines a <c>LIST</c> of <paramref name="directory"/> returns.</summary>
        Task<IReadOnlyList<string>> ListAsync(string directory, CancellationToken cancellationToken);

        /// <summary>Sends a command that opens no data connection (<c>CWD</c>, <c>SIZE</c>, <c>MDTM</c>) and returns its reply, negative or not.</summary>
        Task<MassiveFtpReply> ExecuteAsync(string command, CancellationToken cancellationToken);

        /// <summary>Opens a binary (<c>TYPE I</c>) read of <paramref name="path"/> from byte <paramref name="offset"/> (<c>REST</c>).</summary>
        Task<Stream> OpenReadAsync(string path, long offset, CancellationToken cancellationToken);
    }

    /// <summary>An FTP reply: its three-digit code and its text.</summary>
    internal readonly record struct MassiveFtpReply(int Code, string Message)
    {
        /// <summary>Whether the reply is a positive completion (2xx).</summary>
        public bool IsPositive => Code is >= 200 and < 300;
    }

    /// <summary>The FTP server answered with a negative reply: 4xx (transient by FTP's own definition) or 5xx (permanent).</summary>
    internal sealed class MassiveFtpReplyException(int code, string message, Exception innerException = null)
        : Exception(message, innerException)
    {
        /// <summary>The reply's code, e.g. 450 or 550.</summary>
        public int Code { get; } = code;

        /// <summary>Whether FTP itself calls the reply transient: a 4xx "transient negative completion".</summary>
        public bool IsTransient => Code is >= 400 and < 500;
    }

    /// <summary>
    /// A transport failure talking to the FTP server. An <see cref="IOException"/>, so a failed READ of a file is
    /// reported by <see cref="DownloadEngine"/>'s stream copy exactly as a dropped HTTP body is; a local disk error
    /// is a plain <see cref="IOException"/> and is never mistaken for one of these.
    /// </summary>
    internal sealed class MassiveFtpTransportException(string message, Exception innerException)
        : IOException(message, innerException);

    /// <summary>The <see cref="IMassiveFtpConnection"/> over FluentFTP: explicit FTPS, anonymous, certificate validation on.</summary>
    internal sealed class FluentFtpConnection : IMassiveFtpConnection
    {
        private readonly AsyncFtpClient _client;

        private FluentFtpConnection(AsyncFtpClient client) => _client = client;

        /// <summary>
        /// Connects to <paramref name="host"/> with explicit TLS (<c>AUTH TLS</c>) and logs in anonymously.
        /// </summary>
        /// <param name="host">The FTP host, e.g. <see cref="MassiveClient.FtpHost"/>.</param>
        /// <param name="port">The control port: 21, or a loopback test server's.</param>
        /// <param name="controlTimeout">How long a command's reply may take.</param>
        /// <param name="dataReadTimeout">How long a read of a file may be silent before FluentFTP gives up.</param>
        /// <param name="cancellationToken">Cancels the connection attempt.</param>
        internal static async Task<IMassiveFtpConnection> ConnectAsync(string host, int port, TimeSpan controlTimeout,
            TimeSpan dataReadTimeout, CancellationToken cancellationToken)
        {
            AsyncFtpClient client = CreateClient(host, port, controlTimeout, dataReadTimeout);
            try
            {
                await client.Connect(cancellationToken).ConfigureAwait(false);
                return new FluentFtpConnection(client);
            }
            catch (Exception e) when (Translate(e) is Exception translated)
            {
                client.Dispose();
                throw translated;
            }
            catch
            {
                client.Dispose();
                throw;
            }
        }

        /// <summary>The client, before it connects. Internal so a test can pin its settings.</summary>
        internal static AsyncFtpClient CreateClient(string host, int port, TimeSpan controlTimeout, TimeSpan dataReadTimeout)
        {
            var client = new AsyncFtpClient(host, "anonymous", "anonymous@", port, CreateConfig(controlTimeout, dataReadTimeout));
            client.ValidateCertificate += (_, e) => e.Accept = AcceptCertificate(e.PolicyErrors);
            return client;
        }

        /// <summary>
        /// The connection's settings: explicit TLS on both channels, every certificate validated, binary transfers,
        /// passive mode, and no retry of FluentFTP's own (the client's retry loop is the only one).
        /// </summary>
        internal static FtpConfig CreateConfig(TimeSpan controlTimeout, TimeSpan dataReadTimeout) => new()
        {
            EncryptionMode = FtpEncryptionMode.Explicit,
            DataConnectionEncryption = true,
            // Never true: FluentFTP's own examples often accept any certificate. The handler accepts only a clean one.
            ValidateAnyCertificate = false,
            DataConnectionType = FtpDataConnectionType.AutoPassive,
            DownloadDataType = FtpDataType.Binary,
            ListingDataType = FtpDataType.Binary,
            ConnectTimeout = Milliseconds(controlTimeout),
            ReadTimeout = Milliseconds(controlTimeout),
            DataConnectionConnectTimeout = Milliseconds(controlTimeout),
            DataConnectionReadTimeout = Milliseconds(dataReadTimeout),
            RetryAttempts = 1,
        };

        /// <summary>Accepts a server certificate only when the platform's own validation found nothing wrong with it.</summary>
        internal static bool AcceptCertificate(SslPolicyErrors errors) => errors == SslPolicyErrors.None;

        private static int Milliseconds(TimeSpan span) =>
            span.TotalMilliseconds >= int.MaxValue ? int.MaxValue : Math.Max(1, (int)span.TotalMilliseconds);

        /// <summary>
        /// Maps a FluentFTP or socket failure onto the seam's three kinds, or returns null for an exception that is
        /// not a failure of the connection (the caller's cancellation, or a bug), which then propagates as it is.
        /// </summary>
        internal static Exception Translate(Exception exception) => exception switch
        {
            OperationCanceledException => null,
            MassiveFtpReplyException or MassiveFtpTransportException or AuthenticationException => exception,
            FtpInvalidCertificateException or FtpSecurityNotAvailableException =>
                new AuthenticationException($"The FTP server's TLS was refused: {exception.Message}", exception),
            FtpCommandException command when int.TryParse(command.CompletionCode, NumberStyles.None, CultureInfo.InvariantCulture, out int code) =>
                new MassiveFtpReplyException(code, command.Message, command),
            FtpException or TimeoutException or IOException or SocketException or ObjectDisposedException =>
                new MassiveFtpTransportException($"The FTP connection failed: {exception.Message}", exception),
            AggregateException aggregate when aggregate.InnerException is Exception inner =>
                Translate(inner) ?? new MassiveFtpTransportException($"The FTP connection failed: {aggregate.Message}", aggregate),
            _ => null,
        };

        public async Task<IReadOnlyList<string>> ListAsync(string directory, CancellationToken cancellationToken)
        {
            try
            {
                return await _client.ExecuteDownloadText("LIST " + directory, cancellationToken).ConfigureAwait(false);
            }
            catch (Exception e) when (Translate(e) is Exception translated)
            {
                throw translated;
            }
        }

        public async Task<MassiveFtpReply> ExecuteAsync(string command, CancellationToken cancellationToken)
        {
            try
            {
                FtpReply reply = await _client.Execute(command, cancellationToken).ConfigureAwait(false);
                if (!int.TryParse(reply.Code, NumberStyles.None, CultureInfo.InvariantCulture, out int code))
                    throw new MassiveFtpTransportException($"The FTP server sent a reply with no code to '{command.Split(' ')[0]}'.", null);
                return new MassiveFtpReply(code, reply.Message ?? string.Empty);
            }
            catch (Exception e) when (Translate(e) is Exception translated)
            {
                throw translated;
            }
        }

        public async Task<Stream> OpenReadAsync(string path, long offset, CancellationToken cancellationToken)
        {
            try
            {
                Stream stream = await _client.OpenRead(path, FtpDataType.Binary, offset, false, cancellationToken).ConfigureAwait(false);
                return new MassiveFtpReadStream(stream);
            }
            catch (Exception e) when (Translate(e) is Exception translated)
            {
                throw translated;
            }
        }

        public void Dispose() => _client.Dispose();
    }

    /// <summary>
    /// A file being read over FTP: every failure of a read becomes a <see cref="MassiveFtpTransportException"/>, so the
    /// shared stream copy reports it as the transport failure it is, and a cancellation stays one. Closing it reads the
    /// transfer's final reply; a failure there is ignored, because the client checks the bytes against the server's
    /// SIZE itself.
    /// </summary>
    internal sealed class MassiveFtpReadStream(Stream inner) : Stream
    {
        public override bool CanRead => true;
        public override bool CanSeek => false;
        public override bool CanWrite => false;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => throw new NotSupportedException(); set => throw new NotSupportedException(); }

        public override int Read(byte[] buffer, int offset, int count)
        {
            try { return inner.Read(buffer, offset, count); }
            catch (Exception e) when (FluentFtpConnection.Translate(e) is Exception translated) { throw AsTransport(translated); }
        }

        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            try { return await inner.ReadAsync(buffer, cancellationToken).ConfigureAwait(false); }
            catch (Exception e) when (FluentFtpConnection.Translate(e) is Exception translated) { throw AsTransport(translated); }
        }

        public override Task<int> ReadAsync(byte[] buffer, int offset, int count, CancellationToken cancellationToken) =>
            ReadAsync(buffer.AsMemory(offset, count), cancellationToken).AsTask();

        private static Exception AsTransport(Exception translated) =>
            translated as MassiveFtpTransportException
            ?? new MassiveFtpTransportException($"The FTP data connection failed: {translated.Message}", translated);

        public override void Flush() { }
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
        public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();

        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                try { inner.Dispose(); }
                catch (Exception e) when (FluentFtpConnection.Translate(e) != null) { }
            }
            base.Dispose(disposing);
        }

        public override async ValueTask DisposeAsync()
        {
            try { await inner.DisposeAsync().ConfigureAwait(false); }
            catch (Exception e) when (FluentFtpConnection.Translate(e) != null) { }
            await base.DisposeAsync().ConfigureAwait(false);
        }
    }

    /// <summary>What one line of an FTP <c>LIST</c> names.</summary>
    internal enum MassiveFtpEntryKind
    {
        File,
        Directory,
        Link,
    }

    /// <summary>One entry of an FTP <c>LIST</c>: its name, kind, size and listed date.</summary>
    internal readonly record struct MassiveFtpEntry(string Name, MassiveFtpEntryKind Kind, long Size, DateTime? Modified);

    /// <summary>
    /// Reads the Unix <c>ls -l</c> lines MassIVE's FTP server (SFTPGo 2.7.1) sends for <c>LIST</c>.
    /// </summary>
    /// <remarks>
    /// Read here rather than through FluentFTP's own listing, because neither of its listings tells the truth about
    /// a link (measured 2026-11-10 on MSV000090383): <c>MLSD</c> reports <c>ccms_peak</c> as a 40-byte FILE, and
    /// FluentFTP's <c>LIST</c> parser drops the line altogether, since SFTPGo marks a link with an upper-case
    /// <c>L</c>. A line this parser cannot read is refused, not skipped, so nothing can vanish from a listing silently.
    /// </remarks>
    internal static class MassiveFtpListing
    {
        /// <summary>A line FluentFTP's <c>ExecuteDownloadText</c> puts before the server's own lines.</summary>
        internal const string FluentFtpCaptureHeader = "Lines captured:";

        private static readonly Regex UnixLine = new(
            @"^(?<type>[-dlL])[-rwxsStT]{9}[+@.]?\s+\d+\s+\S+\s+\S+\s+(?<size>\d+)\s+(?<month>[A-Za-z]{3})\s+(?<day>\d{1,2})\s+(?<yearOrTime>\d{4}|\d{1,2}:\d{2})\s(?<name>.+)$",
            RegexOptions.CultureInvariant);

        private static readonly Regex TotalLine = new(@"^total\s+\d+$", RegexOptions.CultureInvariant);

        /// <summary>
        /// The entries of one directory's <c>LIST</c>.
        /// </summary>
        /// <param name="lines">The raw lines.</param>
        /// <param name="directory">The directory listed, for the message.</param>
        /// <param name="now">The time a date with no year is read against (<c>ls</c> shows the time, not the year, for an entry from the last six months).</param>
        /// <exception cref="MzLibException">A line is not one this parser can read.</exception>
        internal static List<MassiveFtpEntry> Parse(IEnumerable<string> lines, string directory, DateTime now)
        {
            var entries = new List<MassiveFtpEntry>();
            foreach (string raw in lines)
            {
                string line = raw?.TrimEnd('\r', '\n');
                if (string.IsNullOrWhiteSpace(line) || line == FluentFtpCaptureHeader || TotalLine.IsMatch(line))
                    continue;

                Match match = UnixLine.Match(line);
                if (!match.Success)
                    throw new MzLibException($"MassIVE's FTP server listed a line in '{directory}' that cannot be read, so the listing would be incomplete: '{line}'.");

                MassiveFtpEntryKind kind = match.Groups["type"].Value switch
                {
                    "d" => MassiveFtpEntryKind.Directory,
                    "l" or "L" => MassiveFtpEntryKind.Link,
                    _ => MassiveFtpEntryKind.File,
                };

                string name = match.Groups["name"].Value;
                // ls -l writes a link as "name -> target"; SFTPGo writes the name alone.
                int arrow = kind == MassiveFtpEntryKind.Link ? name.IndexOf(" -> ", StringComparison.Ordinal) : -1;
                if (arrow > 0)
                    name = name.Substring(0, arrow);

                entries.Add(new MassiveFtpEntry(name, kind,
                    long.Parse(match.Groups["size"].Value, NumberStyles.None, CultureInfo.InvariantCulture),
                    ListedDate(match.Groups["month"].Value, match.Groups["day"].Value, match.Groups["yearOrTime"].Value, now)));
            }
            return entries;
        }

        /// <summary>
        /// The date an <c>ls -l</c> line shows: "Sep 24  2024" (day only) or "Nov 10 12:30" (within the last six
        /// months, in the year that keeps it from lying in the future).
        /// </summary>
        private static DateTime? ListedDate(string month, string day, string yearOrTime, DateTime now)
        {
            if (!DateTime.TryParseExact($"{month} {day} 2000", "MMM d yyyy", CultureInfo.InvariantCulture, DateTimeStyles.None, out DateTime monthDay))
                return null;

            if (yearOrTime.Contains(':'))
            {
                if (!TimeSpan.TryParseExact(yearOrTime, @"h\:mm", CultureInfo.InvariantCulture, out TimeSpan time))
                    return null;
                int year = now.Year;
                if (!TryDate(year, monthDay, out DateTime date) || date.Add(time) > now.AddDays(1))
                {
                    year--;
                    if (!TryDate(year, monthDay, out date))
                        return null;
                }
                return date.Add(time);
            }

            return TryDate(int.Parse(yearOrTime, CultureInfo.InvariantCulture), monthDay, out DateTime day0) ? day0 : null;
        }

        private static bool TryDate(int year, DateTime monthDay, out DateTime date)
        {
            date = default;
            if (monthDay.Month == 2 && monthDay.Day == 29 && !DateTime.IsLeapYear(year))
                return false;
            date = new DateTime(year, monthDay.Month, monthDay.Day, 0, 0, 0, DateTimeKind.Unspecified);
            return true;
        }
    }
}
