using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Security.Authentication;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using MzLibUtil;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A client for MassIVE (https://massive.ucsd.edu) that makes MassIVE proteomics datasets available to mzLib. Given
    /// a MassIVE accession ("MSV000096226") or a ProteomeXchange one MassIVE hosts ("PXD057243"), it returns the
    /// dataset's record, the complete listing of its files with their sizes, and the files' bytes.
    /// </summary>
    /// <remarks>
    /// The record comes from MassIVE's PROXI service over HTTPS. The listing and the files come from MassIVE's FTP
    /// server, <see cref="FtpHost"/>, over explicit FTPS with certificate validation on (decisions D1-D6, 2026-10-18):
    /// MassIVE's HTTPS file endpoint rewrites text files, cannot resume, and rate-limits after a few requests.
    /// <para>
    /// The failure contract is PRIDE's (<see cref="PrideArchiveClient"/>), carried over to FTP:
    /// <list type="bullet">
    /// <item>an outage -- no answer, a dropped or stalled connection, a 5xx from PROXI, a transient (4xx) FTP reply --
    /// is an <see cref="HttpRequestException"/>, even over FTP, so "the service is down" is one type for every
    /// repository client. It is retried first (3 times, after 5, 20 and 60 s), and once the retries are spent,
    /// <see cref="ServiceFailure.IsTransient"/> says so;</item>
    /// <item>MassIVE answering, but with something the contract forbids -- no such dataset, a record with no MassIVE
    /// identifier, a permanent (5xx) FTP reply, a file the listing names but the server no longer has, a certificate
    /// that fails validation -- is an <see cref="MzLibException"/>, never retried.</item>
    /// </list>
    /// An FTP reply code is never read as an HTTP status: FTP's 4xx is transient and its 5xx permanent, the
    /// opposite of HTTP's split. Messages name a file by its path below the dataset and "MassIVE's FTP server",
    /// never by a URL or a login.
    /// </para>
    /// <para>
    /// Holds one <see cref="HttpClient"/> and one FTPS connection, opened on first use and reused; files are fetched
    /// one at a time (one transfer per host, PXReprise 019). Use one instance per unit of work and dispose it.
    /// </para>
    /// </remarks>
    public sealed class MassiveClient : IDisposable
    {
        /// <summary>The base address of MassIVE's PROXI service.</summary>
        public const string DefaultProxiBaseAddress = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/";

        /// <summary>
        /// The base address of ProteomeCentral's PROXI service, asked only to name the repository that hosts a
        /// ProteomeXchange accession MassIVE does not know. Internal: ProteomeCentral is a service of its own, and its
        /// address is not this client's to publish.
        /// </summary>
        internal const string ProteomeCentralProxiBaseAddress = "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/";

        /// <summary>
        /// MassIVE's FTP server. Every record names <c>ftp://massive.ucsd.edu/...</c> as its FTP location, but that host
        /// does not answer FTP (2026-11-10, twice); this one does, and its certificate names it. Only the path is taken
        /// from the record.
        /// </summary>
        public const string FtpHost = "massive-ftp.ucsd.edu";

        /// <summary>The FTP control port: explicit FTPS starts in the clear on 21 and upgrades with <c>AUTH TLS</c>.</summary>
        internal const int FtpPort = 21;

        // MSV000096226, PXD057243, RMSV000000001: letters, then digits. Anything else is refused before any request.
        private static readonly Regex AccessionPattern = new(@"^[A-Za-z]{2,5}\d{5,10}$", RegexOptions.CultureInvariant);

        // One safe folder name of a record's FTP location ("v07", "MSV000096226").
        private static readonly Regex FtpLocationSegment = new(@"^[A-Za-z0-9_.-]+$", RegexOptions.CultureInvariant);

        /// <summary>A hard cap on how deep the listing descends; deeper means a broken listing, not a real dataset.</summary>
        private const int MaxDirectoryDepth = 64;

        private readonly HttpClient _httpClient;
        private readonly bool _ownsHttpClient;
        private bool _disposed;

        // PROXI over HTTPS, and the retry loop and stream copy every FTP operation goes through.
        private readonly DownloadEngine _engine;

        // Asked once, never retried: it only names the host in a message.
        private readonly DownloadEngine _proteomeCentral;

        private readonly Func<CancellationToken, Task<IMassiveFtpConnection>> _connectFtp;
        private readonly SemaphoreSlim _ftpGate = new(1, 1);
        private IMassiveFtpConnection _ftp;

        /// <summary>
        /// How long a download waits for the NEXT bytes of a file before abandoning the transfer, as
        /// <see cref="PrideArchiveClient.BodyStallTimeout"/> does (two minutes). Each read gets a fresh window, so this
        /// bounds silence, not duration. A stall is an <see cref="HttpRequestException"/>, and retried.
        /// </summary>
        public TimeSpan BodyStallTimeout { get => _engine.BodyStallTimeout; init => _engine.BodyStallTimeout = value; }

        /// <summary>How long an FTP command's reply, or a connection, may take. 100 s, as PRIDE's HTTP timeout.</summary>
        internal TimeSpan FtpCommandTimeout { get; init; } = TimeSpan.FromSeconds(100);

        /// <summary>How many times a transient failure is retried after the first attempt (3). Tests set 0 to pin a single attempt.</summary>
        internal int MaxRetries { get => _engine.MaxRetries; init => _engine.MaxRetries = value; }

        /// <summary>Waits out a backoff. Tests replace it to record the schedule instead of sleeping through it.</summary>
        internal Func<TimeSpan, CancellationToken, Task> RetryDelay { get => _engine.RetryDelay; init => _engine.RetryDelay = value; }

        /// <summary>The time a listed date with no year is read against. Tests fix it.</summary>
        internal Func<DateTime> Clock { get; init; } = () => DateTime.UtcNow;

        /// <summary>Creates a client with its own <see cref="HttpClient"/>.</summary>
        public MassiveClient()
            : this(new HttpClient { Timeout = TimeSpan.FromSeconds(100) }, ownsHttpClient: true, connectFtp: null)
        {
        }

        /// <summary>
        /// Creates a client over a supplied <see cref="HttpClient"/>, used for PROXI. The caller retains ownership: the
        /// supplied client is NOT disposed by <see cref="Dispose"/>. Requests are absolute, so its base address is unused.
        /// </summary>
        /// <param name="httpClient">The HTTP client to use. Must not be null.</param>
        public MassiveClient(HttpClient httpClient)
            : this(httpClient, ownsHttpClient: false, connectFtp: null)
        {
        }

        /// <summary>Creates a client whose FTP connections come from <paramref name="connectFtp"/>: the offline tests' seam.</summary>
        internal MassiveClient(HttpClient httpClient, Func<CancellationToken, Task<IMassiveFtpConnection>> connectFtp)
            : this(httpClient, ownsHttpClient: false, connectFtp ?? throw new ArgumentNullException(nameof(connectFtp)))
        {
        }

        private MassiveClient(HttpClient httpClient, bool ownsHttpClient, Func<CancellationToken, Task<IMassiveFtpConnection>> connectFtp)
        {
            _httpClient = httpClient ?? throw new ArgumentNullException(nameof(httpClient));
            _ownsHttpClient = ownsHttpClient;
            _engine = new DownloadEngine(_httpClient, "MassIVE");
            _proteomeCentral = new DownloadEngine(_httpClient, "ProteomeCentral") { MaxRetries = 0 };
            // Read when the connection is made, not now: an object initializer runs after this constructor. FluentFTP's
            // own read deadline sits above the stall window, so a stalled file is reported by the stall window.
            _connectFtp = connectFtp ?? (ct => FluentFtpConnection.ConnectAsync(FtpHost, FtpPort, FtpCommandTimeout,
                BodyStallTimeout + TimeSpan.FromSeconds(30), ct));
        }

        // ---- the record ---------------------------------------------------------------------------------------

        /// <summary>
        /// Returns a MassIVE dataset's record. Use <see cref="TryGetDatasetAsync"/> instead when the accession may
        /// simply not exist.
        /// </summary>
        /// <remarks>
        /// MassIVE answers a ProteomeXchange accession hosted elsewhere with the same 404, byte for byte, as one that
        /// does not exist (PXD000001 and PXD999999, 2026-11-10). So for a PXD MassIVE does not know, ProteomeCentral is
        /// asked once which repository hosts it, and the message names it ("hosted at PRIDE, not MassIVE"). If
        /// ProteomeCentral cannot be asked, the message says that instead; the answer is MassIVE's either way.
        /// </remarks>
        /// <param name="accession">A MassIVE accession ("MSV000096226") or a ProteomeXchange one ("PXD057243").</param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns>The dataset's record. Never null.</returns>
        /// <exception cref="ArgumentException">The accession is blank, or is not letters followed by digits.</exception>
        /// <exception cref="MzLibException">MassIVE has no record for the accession (the message names the repository that hosts a PXD, when ProteomeCentral says), or its record is not one.</exception>
        /// <exception cref="HttpRequestException">MassIVE's PROXI was unreachable or answered a status other than 404, after the retries.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<MassiveDataset> GetDatasetAsync(string accession, CancellationToken cancellationToken = default)
        {
            (bool found, MassiveDataset dataset) = await TryGetDatasetAsync(accession, cancellationToken).ConfigureAwait(false);
            if (found)
                return dataset;

            if (accession.StartsWith("PXD", StringComparison.OrdinalIgnoreCase))
                throw new MzLibException(await DescribeMissingPxdAsync(accession, cancellationToken).ConfigureAwait(false));
            throw new MzLibException($"MassIVE has no dataset with accession '{accession}'.");
        }

        /// <summary>
        /// Attempts to fetch a MassIVE dataset's record, reporting an accession MassIVE does not know as a value rather
        /// than an exception.
        /// </summary>
        /// <remarks>
        /// <c>Found</c> is false for one reason only: MassIVE's PROXI answered 404. That includes a ProteomeXchange
        /// accession another repository hosts; <see cref="GetDatasetAsync"/> names the host, this does not ask.
        /// </remarks>
        /// <param name="accession">A MassIVE accession ("MSV000096226") or a ProteomeXchange one ("PXD057243").</param>
        /// <param name="cancellationToken">Cancels the fetch.</param>
        /// <returns><c>Found</c> and the record, or <c>false</c> and null. The record is never null when <c>Found</c> is true.</returns>
        /// <exception cref="ArgumentException">The accession is blank, or is not letters followed by digits.</exception>
        /// <exception cref="MzLibException">MassIVE answered successfully with something that is not a dataset record.</exception>
        /// <exception cref="HttpRequestException">MassIVE's PROXI was unreachable or answered a status other than 404, after the retries.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<(bool Found, MassiveDataset Dataset)> TryGetDatasetAsync(string accession,
            CancellationToken cancellationToken = default)
        {
            ValidateAccession(accession);
            cancellationToken.ThrowIfCancellationRequested();

            string requestUri = $"{DefaultProxiBaseAddress}datasets/{Uri.EscapeDataString(accession)}";
            string described = $"MassIVE's PROXI record for '{accession}'";
            using HttpResponseMessage response = await _engine.WithRetryAsync(async () =>
            {
                HttpResponseMessage attempt = await _engine.GetAsync(requestUri, HttpCompletionOption.ResponseContentRead,
                    cancellationToken, described).ConfigureAwait(false);

                // 404 is the one expected "no"; every other failure still throws.
                if (attempt.IsSuccessStatusCode || attempt.StatusCode == HttpStatusCode.NotFound)
                    return attempt;

                using (attempt)
                    throw DownloadEngine.StatusFailure(attempt,
                        $"MassIVE PROXI request failed with status {(int)attempt.StatusCode} {attempt.ReasonPhrase} for '{accession}'.");
            }, _engine.HostOf(requestUri), cancellationToken).ConfigureAwait(false);

            if (response.StatusCode == HttpStatusCode.NotFound)
                return (false, null);

            // The bytes, not ReadAsStringAsync: that follows the ISO-8859-1 label MassIVE puts on every record, and
            // half of them are UTF-8.
            byte[] body = await response.Content.ReadAsByteArrayAsync(cancellationToken).ConfigureAwait(false);
            if (body.Length == 0)
                throw new MzLibException($"MassIVE's PROXI returned an empty body for '{accession}'.");
            return (true, MassiveDataset.FromProxi(body, accession));
        }

        /// <summary>The message for a PXD MassIVE does not know, naming its host when ProteomeCentral says.</summary>
        private async Task<string> DescribeMissingPxdAsync(string accession, CancellationToken cancellationToken)
        {
            string requestUri = $"{ProteomeCentralProxiBaseAddress}datasets/{Uri.EscapeDataString(accession)}";
            try
            {
                using HttpResponseMessage response = await _proteomeCentral.GetAsync(requestUri,
                    HttpCompletionOption.ResponseContentRead, cancellationToken, $"ProteomeCentral's record for '{accession}'").ConfigureAwait(false);
                if (response.StatusCode == HttpStatusCode.NotFound)
                    return $"Neither MassIVE nor ProteomeCentral has a dataset with accession '{accession}'.";
                if (response.IsSuccessStatusCode)
                {
                    string json = await response.Content.ReadAsStringAsync(cancellationToken).ConfigureAwait(false);
                    JToken hostToken = (JToken.Parse(json) as JObject)?["datasetSummary"]?["hostingRepository"];
                    string host = hostToken?.Type == JTokenType.String ? (string)hostToken : null;
                    if (!string.IsNullOrWhiteSpace(host))
                        return string.Equals(host, "MassIVE", StringComparison.OrdinalIgnoreCase)
                            ? $"ProteomeCentral says MassIVE hosts '{accession}', but MassIVE's PROXI has no record for it."
                            : $"'{accession}' is hosted at {host}, not MassIVE (ProteomeCentral's record says so), so MassIVE has no dataset for it.";
                }
            }
            catch (HttpRequestException) { }
            catch (JsonException) { }

            return $"MassIVE has no dataset with accession '{accession}', and ProteomeCentral could not be asked where it is hosted.";
        }

        private static void ValidateAccession(string accession)
        {
            if (string.IsNullOrWhiteSpace(accession))
                throw new ArgumentException("A MassIVE or ProteomeXchange accession is required.", nameof(accession));
            if (!AccessionPattern.IsMatch(accession))
                throw new ArgumentException(
                    $"'{accession}' is not a MassIVE or ProteomeXchange accession (e.g. \"MSV000096226\" or \"PXD057243\").", nameof(accession));
        }

        // ---- the listing --------------------------------------------------------------------------------------

        /// <summary>
        /// Lists every entry below a dataset's folder on MassIVE's FTP server, subdirectories and updates included, with
        /// each file's exact size, before anything is downloaded.
        /// </summary>
        /// <remarks>
        /// <b>Never empty.</b> A dataset MassIVE has a record for always has files, so a root that lists nothing is
        /// MassIVE's FTP server failing to list, not an empty dataset: it is retried as a transport failure, as PRIDE's
        /// empty FTP root is (PXD058248), and then thrown.
        /// <para>
        /// <b>Links</b> are listed with <see cref="MassiveFile.IsLink"/> set and are not followed: MassIVE's FTP server
        /// will not open them. <c>ccms_peak</c> is one on every dataset checked (2026-11-10), so MassIVE's own converted
        /// copies are never in the listing.
        /// </para>
        /// </remarks>
        /// <param name="accession">A MassIVE accession, or a ProteomeXchange one MassIVE hosts.</param>
        /// <param name="cancellationToken">Cancels the (multi-request) listing.</param>
        /// <returns>Every file and link below the dataset's folder, in listing order. Never null, never empty.</returns>
        /// <exception cref="ArgumentException">The accession is blank, or is not letters followed by digits.</exception>
        /// <exception cref="MzLibException">As for <see cref="GetDatasetAsync"/>; or the record names no usable FTP location, the server has no such folder, a listing line cannot be read, the listing nests implausibly deep, the dataset lists only empty folders, or a permanent FTP reply refused a listing.</exception>
        /// <exception cref="HttpRequestException">PROXI or the FTP server failed transiently on every attempt, the root listed nothing on every attempt included.</exception>
        /// <exception cref="OperationCanceledException">The operation was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<IReadOnlyList<MassiveFile>> GetDatasetFilesAsync(string accession, CancellationToken cancellationToken = default)
        {
            MassiveDataset dataset = await GetDatasetAsync(accession, cancellationToken).ConfigureAwait(false);
            string root = DatasetRootOf(dataset);

            var files = new List<MassiveFile>();
            await CollectAsync(dataset.MassiveAccession, root, relative: "", depth: 0, files, cancellationToken).ConfigureAwait(false);
            if (files.Count == 0)
                throw new MzLibException($"MassIVE's FTP server lists only empty folders below '{root}' for {dataset.MassiveAccession}.");
            return files;
        }

        /// <summary>The dataset's folder on the FTP server, "/" + its record's FTP location, checked to be a plain path.</summary>
        private static string DatasetRootOf(MassiveDataset dataset)
        {
            string location = dataset.FtpLocation;
            string[] segments = location?.Split('/') ?? Array.Empty<string>();
            if (segments.Length == 0 || segments.Any(s => s is "." or ".." || !FtpLocationSegment.IsMatch(s)))
                throw new MzLibException(location == null
                    ? $"MassIVE's PROXI record for {dataset.MassiveAccession} names no FTP location, so its files cannot be listed."
                    : $"MassIVE's PROXI record for {dataset.MassiveAccession} names an FTP location that is not a plain folder path: '{location}'.");
            return "/" + string.Join("/", segments);
        }

        /// <summary>Lists one folder, adds its files and links to <paramref name="files"/>, and descends into its folders.</summary>
        private async Task CollectAsync(string accession, string root, string relative, int depth, List<MassiveFile> files,
            CancellationToken cancellationToken)
        {
            cancellationToken.ThrowIfCancellationRequested();
            if (depth > MaxDirectoryDepth)
                throw new MzLibException(
                    $"MassIVE's FTP listing of {accession} nests deeper than {MaxDirectoryDepth} folders at '{relative}'; the listing may be cyclic.");

            string directory = relative.Length == 0 ? root : root + "/" + relative;
            string described = relative.Length == 0 ? $"the listing of {accession}'s folder" : $"the listing of '{relative}' in {accession}";

            List<MassiveFtpEntry> entries = await _engine.WithRetryAsync(() => UseFtpAsync(async ftp =>
            {
                if (depth == 0)
                {
                    // MassIVE's FTP server answers a LIST of a missing folder with 450, FTP's "try again later", so a
                    // folder the record names but the server lacks would be retried for 85 s. CWD says 550 for it.
                    MassiveFtpReply cwd = await ftp.ExecuteAsync("CWD " + directory, cancellationToken).ConfigureAwait(false);
                    if (cwd.Code == 550)
                        throw new MzLibException(
                            $"MassIVE's FTP server has no folder '{directory}', which {accession}'s PROXI record names as its FTP location: {cwd.Message}");
                    ThrowUnlessPositive(cwd, described);
                }

                IReadOnlyList<string> lines = await ftp.ListAsync(directory, cancellationToken).ConfigureAwait(false);
                List<MassiveFtpEntry> listed = MassiveFtpListing.Parse(lines, directory, Clock());
                if (depth == 0 && listed.Count == 0)
                    throw new HttpRequestException(
                        $"MassIVE's FTP listing of {accession}'s folder came back empty although MassIVE has a record for it.");
                return listed;
            }, described, cancellationToken), FtpHost, cancellationToken).ConfigureAwait(false);

            foreach (MassiveFtpEntry entry in entries)
            {
                if (entry.Name.Length == 0 || entry.Name is "." or ".." || entry.Name.Contains('/'))
                    continue;

                string path = relative.Length == 0 ? entry.Name : relative + "/" + entry.Name;
                switch (entry.Kind)
                {
                    case MassiveFtpEntryKind.Directory:
                        await CollectAsync(accession, root, path, depth + 1, files, cancellationToken).ConfigureAwait(false);
                        break;
                    case MassiveFtpEntryKind.Link:
                        files.Add(new MassiveFile(accession, root, path, sizeBytes: null, entry.Modified, isLink: true));
                        break;
                    default:
                        files.Add(new MassiveFile(accession, root, path, entry.Size, entry.Modified));
                        break;
                }
            }
        }

        // ---- downloads ----------------------------------------------------------------------------------------

        /// <summary>
        /// Downloads one file from MassIVE's FTP server below <paramref name="destinationDirectory"/>, keeping its path
        /// below the dataset ("raw/a/x.raw" lands in "&lt;destination&gt;/raw/a/x.raw"), and checks its size.
        /// </summary>
        /// <remarks>
        /// <b>Always size-checked.</b> Before the transfer the server's <c>SIZE</c> must equal the listing's. If it does
        /// not, the file's folder is listed again once: a file no longer there, or a size the new listing does not agree
        /// with either, is an <see cref="MzLibException"/>; a size it agrees with is the one checked. After the
        /// transfer the bytes written must equal it, or the attempt is retried from byte 0.
        /// <para>
        /// <b>Binary</b> (<c>TYPE I</c>), so a text file arrives byte for byte; FTP's ASCII mode would rewrite its line
        /// ends, as MassIVE's HTTPS endpoint does (sdrf 028).
        /// </para>
        /// <para>
        /// <b>Resumable.</b> The bytes go to "&lt;file&gt;.partial", beside a ".partial.validator" holding the
        /// server's <c>SIZE</c> and <c>MDTM</c> when the transfer began. A later attempt, in this call or a later one,
        /// resumes with <c>REST</c> only while both still match; otherwise it starts again from byte 0. Cancelling, or a
        /// transient failure, keeps both files for the next call; any other failure removes them.
        /// </para>
        /// </remarks>
        /// <param name="file">A file from <see cref="GetDatasetFilesAsync"/>. Must not be null or a link.</param>
        /// <param name="destinationDirectory">The directory to write below; created if it does not exist.</param>
        /// <param name="overwrite">
        /// When true (the default) an existing destination file is replaced. When false, an existing destination file
        /// of the listed size is left untouched and no request is made; one of another size is downloaded again.
        /// </param>
        /// <param name="cancellationToken">Cancels the download.</param>
        /// <returns>The full path of the written (or already-present) file.</returns>
        /// <exception cref="ArgumentNullException">The file is null.</exception>
        /// <exception cref="ArgumentException">The destination directory is blank, or a segment of the file's path cannot be a local file name (a "..", a "." or an empty segment, or an invalid character), so it could escape the destination.</exception>
        /// <exception cref="NotSupportedException">The entry is a link, which MassIVE's FTP server will not open.</exception>
        /// <exception cref="MzLibException">The server no longer has the file, its size disagrees with the listing after listing it again, a permanent FTP reply refused it, or the TLS certificate failed validation.</exception>
        /// <exception cref="HttpRequestException">The FTP server failed transiently on every attempt: no answer, a dropped or stalled transfer, a 4xx reply, or a file that came up short.</exception>
        /// <exception cref="OperationCanceledException">The download was cancelled via <paramref name="cancellationToken"/>.</exception>
        public async Task<string> DownloadFileAsync(MassiveFile file, string destinationDirectory, bool overwrite = true,
            CancellationToken cancellationToken = default)
        {
            if (file == null)
                throw new ArgumentNullException(nameof(file));
            if (string.IsNullOrWhiteSpace(destinationDirectory))
                throw new ArgumentException("A destination directory is required.", nameof(destinationDirectory));
            if (file.IsLink)
                throw new NotSupportedException(
                    $"'{file.Path}' is a link on MassIVE's FTP server, which will not open links, so it cannot be downloaded.");

            string destinationPath = LocalPathOf(file, destinationDirectory);
            long expectedSize = file.SizeBytes.Value;

            // Cheap resume: a destination file already of the listed size is kept, with no request made.
            if (!overwrite && File.Exists(destinationPath) && new FileInfo(destinationPath).Length == expectedSize)
                return destinationPath;

            string described = $"'{file.Path}' from MassIVE's FTP server";
            Directory.CreateDirectory(Path.GetDirectoryName(destinationPath));
            string partialPath = destinationPath + ".partial";
            string validatorPath = partialPath + PrideArchiveClient.ValidatorSuffix;

            FtpValidator onDisk = FtpValidator.Read(partialPath, validatorPath);
            bool relisted = false;

            try
            {
                await _engine.WithRetryAsync(() => UseFtpAsync(async ftp =>
                {
                    long? sizeNow = await SizeAsync(ftp, file.ServerPath, described, cancellationToken).ConfigureAwait(false);
                    if (sizeNow != expectedSize)
                    {
                        if (relisted)
                            throw new MzLibException(sizeNow == null
                                ? $"MassIVE's FTP server no longer has {described}."
                                : $"{described} changed while it was being downloaded: the server's SIZE was {expectedSize} bytes and is now {sizeNow}.");

                        relisted = true;
                        MassiveFtpEntry? fresh = await RelistAsync(ftp, file, cancellationToken).ConfigureAwait(false);
                        if (fresh == null)
                            throw new MzLibException($"MassIVE's FTP server no longer has {described}: its folder, listed again, does not hold it.");
                        if (sizeNow == null)
                            throw new HttpRequestException(
                                $"MassIVE's FTP server answered SIZE for {described} with 550, yet still lists it.");
                        if (fresh.Value.Size != sizeNow.Value)
                            throw new MzLibException(
                                $"The listing of {described} is stale: it gave {expectedSize} bytes, the server's SIZE gives {sizeNow}, and listing it again gives {fresh.Value.Size}.");
                        expectedSize = sizeNow.Value;
                    }

                    string modified = await ModifiedAsync(ftp, file.ServerPath, described, cancellationToken).ConfigureAwait(false)
                        ?? throw new MzLibException($"MassIVE's FTP server no longer has {described}.");
                    var current = new FtpValidator(expectedSize, modified);

                    long offset = 0;
                    if (current.Equals(onDisk) && File.Exists(partialPath))
                    {
                        long have = new FileInfo(partialPath).Length;
                        if (have <= expectedSize)
                            offset = have;
                    }
                    if (offset == 0)
                    {
                        current.Write(validatorPath);
                        onDisk = current;
                    }

                    if (offset < expectedSize || offset == 0)
                    {
                        using (var fileStream = new FileStream(partialPath, offset > 0 ? FileMode.Append : FileMode.Create,
                                   FileAccess.Write, FileShare.None))
                        using (Stream data = await ftp.OpenReadAsync(file.ServerPath, offset, cancellationToken).ConfigureAwait(false))
                        {
                            await _engine.CopyUntilStalledAsync(data, fileStream, described, inspector: null, cancellationToken).ConfigureAwait(false);
                        }
                    }

                    long actual = new FileInfo(partialPath).Length;
                    if (actual != expectedSize)
                    {
                        // Not trusted for a resume either: the next attempt starts again from byte 0.
                        DownloadEngine.DeleteQuietly(partialPath);
                        DownloadEngine.DeleteQuietly(validatorPath);
                        onDisk = null;
                        throw new HttpRequestException(
                            $"The MassIVE download of {described} came to {actual} bytes where the server's SIZE gives {expectedSize}.");
                    }
                    return true;
                }, described, cancellationToken), FtpHost, cancellationToken).ConfigureAwait(false);

                File.Move(partialPath, destinationPath, overwrite: true);
                DownloadEngine.DeleteQuietly(validatorPath);
            }
            catch (HttpRequestException e) when (_engine.IsTransient(e, FtpHost) && onDisk != null && File.Exists(partialPath))
            {
                // Kept on purpose: the server failed, not the file, and the validator lets the next call resume.
                throw;
            }
            catch (OperationCanceledException) when (onDisk != null && File.Exists(partialPath))
            {
                // Kept on purpose: the caller stopped the transfer, not the file (PXReprise 019).
                throw;
            }
            catch
            {
                DownloadEngine.DeleteQuietly(partialPath);
                DownloadEngine.DeleteQuietly(validatorPath);
                throw;
            }

            return destinationPath;
        }

        /// <summary>
        /// Downloads a dataset's files, as <see cref="DownloadFileAsync"/> downloads each, one at a time in listing
        /// order. Links are never downloaded: they are not files MassIVE's FTP server will serve.
        /// </summary>
        /// <param name="accession">A MassIVE accession, or a ProteomeXchange one MassIVE hosts.</param>
        /// <param name="destinationDirectory">The directory to write below; created if it does not exist.</param>
        /// <param name="filter">An optional predicate choosing which files to download. When null, every file is downloaded.</param>
        /// <param name="overwrite">Passed through to each download.</param>
        /// <param name="cancellationToken">Cancels between and during file downloads.</param>
        /// <returns>The full paths of the downloaded files, in listing order. Empty if the filter chose none.</returns>
        /// <exception cref="ArgumentException">The accession or destination directory is blank, or the accession is not letters followed by digits.</exception>
        /// <exception cref="MzLibException">As for <see cref="GetDatasetFilesAsync"/> and <see cref="DownloadFileAsync"/>.</exception>
        /// <exception cref="HttpRequestException">As for <see cref="GetDatasetFilesAsync"/> and <see cref="DownloadFileAsync"/>.</exception>
        public async Task<IReadOnlyList<string>> DownloadDatasetFilesAsync(string accession, string destinationDirectory,
            Func<MassiveFile, bool> filter = null, bool overwrite = true, CancellationToken cancellationToken = default)
        {
            if (string.IsNullOrWhiteSpace(destinationDirectory))
                throw new ArgumentException("A destination directory is required.", nameof(destinationDirectory));

            IReadOnlyList<MassiveFile> files = await GetDatasetFilesAsync(accession, cancellationToken).ConfigureAwait(false);
            var paths = new List<string>();
            foreach (MassiveFile file in files.Where(f => !f.IsLink && (filter == null || filter(f))))
            {
                cancellationToken.ThrowIfCancellationRequested();
                paths.Add(await DownloadFileAsync(file, destinationDirectory, overwrite, cancellationToken).ConfigureAwait(false));
            }
            return paths;
        }

        /// <summary>
        /// Where <paramref name="file"/> is written below <paramref name="destinationDirectory"/>. Every segment of its
        /// path must be a plain local file name, so a path from the server cannot write outside the destination.
        /// </summary>
        private static string LocalPathOf(MassiveFile file, string destinationDirectory)
        {
            char[] invalid = Path.GetInvalidFileNameChars();
            string[] segments = file.Path.Split('/');
            foreach (string segment in segments)
            {
                if (segment.Length == 0 || segment is "." or ".." || segment.IndexOfAny(invalid) >= 0
                    || segment.Contains('\\') || segment.Contains(':') || Path.IsPathRooted(segment))
                    throw new ArgumentException(
                        $"The MassIVE path '{file.Path}' has a segment that cannot be a local file name; refusing to write outside the destination directory.",
                        nameof(file));
            }

            string root = Path.GetFullPath(destinationDirectory);
            string path = Path.GetFullPath(Path.Combine(new[] { root }.Concat(segments).ToArray()));
            if (!path.StartsWith(Path.TrimEndingDirectorySeparator(root) + Path.DirectorySeparatorChar, StringComparison.Ordinal))
                throw new ArgumentException(
                    $"The MassIVE path '{file.Path}' resolves outside the destination directory.", nameof(file));
            return path;
        }

        /// <summary>The server's SIZE of <paramref name="serverPath"/>, or null when it answers 550 (no such file).</summary>
        private static async Task<long?> SizeAsync(IMassiveFtpConnection ftp, string serverPath, string described, CancellationToken cancellationToken)
        {
            MassiveFtpReply reply = await ftp.ExecuteAsync("SIZE " + serverPath, cancellationToken).ConfigureAwait(false);
            if (reply.Code == 550)
                return null;
            ThrowUnlessPositive(reply, described);
            if (!long.TryParse(reply.Message.Trim(), NumberStyles.None, CultureInfo.InvariantCulture, out long size))
                throw new MzLibException($"MassIVE's FTP server answered SIZE for {described} with '{reply.Message}', which is not a size.");
            return size;
        }

        /// <summary>The server's MDTM stamp of <paramref name="serverPath"/> as sent ("20150212140541"), or null when it answers 550.</summary>
        private static async Task<string> ModifiedAsync(IMassiveFtpConnection ftp, string serverPath, string described, CancellationToken cancellationToken)
        {
            MassiveFtpReply reply = await ftp.ExecuteAsync("MDTM " + serverPath, cancellationToken).ConfigureAwait(false);
            if (reply.Code == 550)
                return null;
            ThrowUnlessPositive(reply, described);
            string stamp = reply.Message.Trim();
            if (stamp.Length == 0 || stamp.Any(char.IsWhiteSpace))
                throw new MzLibException($"MassIVE's FTP server answered MDTM for {described} with '{reply.Message}', which is not a time stamp.");
            return stamp;
        }

        /// <summary>Lists <paramref name="file"/>'s folder again and returns its entry, or null when it is not there as a file.</summary>
        private async Task<MassiveFtpEntry?> RelistAsync(IMassiveFtpConnection ftp, MassiveFile file, CancellationToken cancellationToken)
        {
            string folder = file.ServerPath.Substring(0, file.ServerPath.LastIndexOf('/'));
            IReadOnlyList<string> lines = await ftp.ListAsync(folder, cancellationToken).ConfigureAwait(false);
            foreach (MassiveFtpEntry entry in MassiveFtpListing.Parse(lines, folder, Clock()))
            {
                if (entry.Kind == MassiveFtpEntryKind.File && string.Equals(entry.Name, file.FileName, StringComparison.Ordinal))
                    return entry;
            }
            return null;
        }

        // ---- the FTP connection -------------------------------------------------------------------------------

        /// <summary>
        /// Runs one FTP operation on the client's connection, opening it if needed, one operation at a time, and maps
        /// a failure of the connection onto the client's contract. After any failure the connection is dropped, and the
        /// next operation opens a new one: a failed command or an abandoned transfer leaves FTP's control channel in a
        /// state nothing can vouch for.
        /// </summary>
        private async Task<T> UseFtpAsync<T>(Func<IMassiveFtpConnection, Task<T>> operation, string described,
            CancellationToken cancellationToken)
        {
            ObjectDisposedException.ThrowIf(_disposed, this);
            await _ftpGate.WaitAsync(cancellationToken).ConfigureAwait(false);
            try
            {
                _ftp ??= await _connectFtp(cancellationToken).ConfigureAwait(false);
                return await operation(_ftp).ConfigureAwait(false);
            }
            catch (Exception e)
            {
                DropFtp();
                Exception mapped = e switch
                {
                    MassiveFtpTransportException => new HttpRequestException(
                        $"MassIVE's FTP server connection failed during {described}.", e),
                    AuthenticationException => new MzLibException(
                        $"MassIVE's FTP server failed the TLS certificate check during {described}, so nothing was sent or received over it.", e),
                    MassiveFtpReplyException reply => ReplyFailure(reply.Code, reply.Message, described, reply),
                    _ => null,
                };
                if (mapped == null)
                    throw;
                throw mapped;
            }
            finally
            {
                _ftpGate.Release();
            }
        }

        /// <summary>
        /// The contract's exception for a negative FTP reply: a 4xx, FTP's "transient negative completion", is a
        /// status-less <see cref="HttpRequestException"/>, so it is retried and then marked transient; a 5xx is FTP's
        /// "permanent negative completion" and an <see cref="MzLibException"/>.
        /// </summary>
        private static Exception ReplyFailure(int code, string message, string described, Exception innerException) =>
            code is >= 400 and < 500
                ? new HttpRequestException($"MassIVE's FTP server answered {code}, FTP's 'try again later', to {described}: {message}", innerException)
                : new MzLibException($"MassIVE's FTP server refused {described} with {code}: {message}", innerException);

        private static void ThrowUnlessPositive(MassiveFtpReply reply, string described)
        {
            if (!reply.IsPositive)
                throw ReplyFailure(reply.Code, reply.Message, described, null);
        }

        private void DropFtp()
        {
            IMassiveFtpConnection ftp = _ftp;
            _ftp = null;
            try { ftp?.Dispose(); }
            catch (Exception) { /* closing a failed connection must not replace the failure that closed it */ }
        }

        /// <summary>Closes the FTP connection, and the <see cref="HttpClient"/> when this client created it.</summary>
        public void Dispose()
        {
            if (_disposed)
                return;
            _disposed = true;
            DropFtp();
            if (_ownsHttpClient)
                _httpClient.Dispose();
        }

        /// <summary>
        /// What a <c>.partial</c> was started under on the FTP server: its <c>SIZE</c> and <c>MDTM</c> then. A resume is
        /// made only while both still match. Kept in the ".partial.validator" sidecar PRIDE's downloads use, as one line,
        /// "ftp size=&lt;bytes&gt; mdtm=&lt;stamp&gt;".
        /// </summary>
        private sealed record FtpValidator(long Size, string Modified)
        {
            private static readonly Regex Line = new(@"^ftp size=(?<size>\d+) mdtm=(?<mdtm>\S+)$", RegexOptions.CultureInvariant);

            public override string ToString() => $"ftp size={Size.ToString(CultureInfo.InvariantCulture)} mdtm={Modified}";

            public void Write(string validatorPath) => File.WriteAllText(validatorPath, ToString());

            /// <summary>
            /// The validator a partial was started under, or null. An orphan -- a partial with no readable validator, or a
            /// validator with no partial -- is deleted, because it can never be resumed.
            /// </summary>
            public static FtpValidator Read(string partialPath, string validatorPath)
            {
                FtpValidator validator = null;
                try
                {
                    if (File.Exists(partialPath) && File.Exists(validatorPath))
                    {
                        Match match = Line.Match(File.ReadAllText(validatorPath).Trim());
                        if (match.Success && long.TryParse(match.Groups["size"].Value, NumberStyles.None, CultureInfo.InvariantCulture, out long size))
                            validator = new FtpValidator(size, match.Groups["mdtm"].Value);
                    }
                }
                catch (IOException) { }
                catch (UnauthorizedAccessException) { }

                if (validator == null)
                {
                    DownloadEngine.DeleteQuietly(partialPath);
                    DownloadEngine.DeleteQuietly(validatorPath);
                }
                return validator;
            }
        }
    }
}
