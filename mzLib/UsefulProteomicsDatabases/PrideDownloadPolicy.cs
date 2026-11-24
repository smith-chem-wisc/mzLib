using System;
using System.Threading;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// How <see cref="PrideArchiveClient.DownloadFileWithPolicyAsync"/> gets one PRIDE file: which routes it tries, how
    /// often, how long it waits between tries, and when it stops using Aspera for a while.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Opt-in. mzLib's own downloads make one Aspera attempt
    /// (<see cref="PrideArchiveClient.DownloadFileViaAsperaAsync"/>) or one HTTPS download with its fixed internal retry
    /// (<see cref="PrideArchiveClient.DownloadFileAsync(PrideArchiveFile, string, PrideFileChecksum, bool, bool, CancellationToken)"/>),
    /// and report what went wrong. This is for a caller that downloads PRIDE continuously and wants mzLib to own that
    /// policy as well.
    /// </para>
    /// <para>
    /// The defaults are PXReprise's, which has downloaded thousands of PRIDE files over both routes since 2026-10-08:
    /// eight tries in all, waits growing by 15 s, three Aspera tries before HTTPS takes over, a damaged Aspera copy
    /// sent to HTTPS (Aspera damaged 10.6% of 1,754 copies to HTTPS's 5.2% of 194, PXReprise 019), three Aspera sessions
    /// at once with starts 10 s apart, and Aspera set aside for 30 min after EBI refuses a login or a connection.
    /// </para>
    /// </remarks>
    public sealed class PrideDownloadPolicy
    {
        /// <summary>How Aspera is run; null (the default) downloads over HTTPS only.</summary>
        public PrideAsperaOptions Aspera { get; init; }

        /// <summary>The most tries, over all routes together. Eight by default.</summary>
        public int Attempts { get; init; } = 8;

        /// <summary>
        /// The wait before the second try; the wait before try n+1 is n times this, so the waits grow by it. 15 s by
        /// default: 15 s, 30 s, 45 s and so on.
        /// </summary>
        public TimeSpan Backoff { get; init; } = TimeSpan.FromSeconds(15);

        /// <summary>How many of a file's tries may go to Aspera before the rest go to HTTPS. Three by default.</summary>
        public int AsperaTries { get; init; } = 3;

        /// <summary>
        /// Whether a file whose Aspera copy came back damaged goes over HTTPS for the rest of its tries. True by default:
        /// Aspera damaged twice as many copies as HTTPS did, and more the larger the file.
        /// </summary>
        public bool HttpsAfterDamagedAsperaCopy { get; init; } = true;

        /// <summary>
        /// How long Aspera is set aside, for every download in the process, after a
        /// <see cref="PrideAsperaFailure.LoginRefused"/> or <see cref="PrideAsperaFailure.CannotConnect"/>: those are
        /// failures of the route, and every file would fail the same way. Files go over HTTPS meanwhile; the first
        /// Aspera try after the pause tests the route again. 30 min by default; zero sets nothing aside.
        /// </summary>
        /// <remarks>
        /// EBI refused this route for about three days from 2026-10-30, and retrying it per file cost each file about a
        /// minute for nothing (PXReprise 032, 033).
        /// </remarks>
        public TimeSpan AsperaPauseAfterRouteFailure { get; init; } = TimeSpan.FromMinutes(30);

        /// <summary>The most Aspera sessions running at once in the process, under any policy. Three by default.</summary>
        public int MaxConcurrentAsperaSessions { get; init; } = 3;

        /// <summary>The least time between two Aspera sessions' starts in the process. 10 s by default.</summary>
        public TimeSpan AsperaStartGap { get; init; } = TimeSpan.FromSeconds(10);

        /// <summary>Receives a <see cref="PrideDownloadTry"/> as each try ends, accepted or not; null (the default) receives none.</summary>
        public IProgress<PrideDownloadTry> Tries { get; init; }

        /// <summary>
        /// Picks where each HTTPS try downloads from. It is given the file and the record of the file's previous try in
        /// this call, over either route (null before the first try). It returns an absolute <c>https</c> URL, or null for
        /// the file's own HTTPS location (<see cref="PrideArchiveExtensions.TryGetHttpsDownloadUrl"/>). Null (the default)
        /// uses the file's own location for every try.
        /// </summary>
        /// <remarks>
        /// <para>
        /// For a caller with more than one HTTPS host for PRIDE's files. PXReprise tries EBI's Globus HTTPS endpoint
        /// first; after a Globus try that fails without a damaged copy, it makes one try from <c>ftp.pride.ebi.ac.uk</c>,
        /// then goes back to Globus (PXReprise 040). The mapping onto another host stays the caller's.
        /// </para>
        /// <para>
        /// It is asked only for a file that has an HTTPS location of its own, and never for an Aspera try. The host of
        /// the URL it returns is the try's <see cref="PrideDownloadTry.Host"/>. A partial copy is resumed from another
        /// host only when its MD5 is checked: the validator it was started under is that host's, and means nothing to
        /// another. Anything but an absolute <c>https</c> URL is an <see cref="ArgumentException"/>, thrown at once, and
        /// an exception it throws is thrown at once too.
        /// </para>
        /// </remarks>
        public Func<PrideArchiveFile, PrideDownloadTry, string> HttpsLocation { get; init; }
    }

    /// <summary>The route one try of <see cref="PrideArchiveClient.DownloadFileWithPolicyAsync"/> took.</summary>
    public enum PrideDownloadRoute
    {
        /// <summary>EBI's Aspera server, through <c>ascp</c>.</summary>
        Aspera,

        /// <summary>HTTPS, from the file's HTTPS location or PRIDE's FTP host over HTTPS.</summary>
        Https
    }

    /// <summary>How one try of <see cref="PrideArchiveClient.DownloadFileWithPolicyAsync"/> ended.</summary>
    public enum PrideDownloadOutcome
    {
        /// <summary>The file arrived and passed its check against PRIDE's checksum row.</summary>
        Accepted,

        /// <summary>The file arrived but failed its check, and was discarded.</summary>
        Damaged,

        /// <summary>The try failed before a whole file arrived.</summary>
        Failed
    }

    /// <summary>One try of <see cref="PrideArchiveClient.DownloadFileWithPolicyAsync"/>, as it ended.</summary>
    /// <param name="FileName">The file's name.</param>
    /// <param name="Try">The try's number, from 1.</param>
    /// <param name="Route">The route it took.</param>
    /// <param name="Host">The host it downloaded from.</param>
    /// <param name="FileSizeBytes">The file's size by PRIDE's checksum row.</param>
    /// <param name="StartedUtc">
    /// When the transfer started, in UTC: after any wait for an Aspera session, so that tries which overlapped can be
    /// counted. For a try that failed before its transfer started, when it failed.
    /// </param>
    /// <param name="Elapsed">From the start of the transfer to its end, not counting waits for an Aspera session.</param>
    /// <param name="BytesReceived">
    /// The bytes this try received, whole or not; null when they could not be counted. An accepted copy received its
    /// size, unless it resumed a partial copy an earlier try left (then only the rest), or was found already whole on
    /// disk (then 0). Over Aspera, the bytes in ascp's folder when it ended.
    /// </param>
    /// <param name="Outcome">How it ended.</param>
    /// <param name="Error">Why it did not succeed; null when <see cref="PrideDownloadOutcome.Accepted"/>.</param>
    public sealed record PrideDownloadTry(string FileName, int Try, PrideDownloadRoute Route, string Host,
        long FileSizeBytes, DateTime StartedUtc, TimeSpan Elapsed, long? BytesReceived, PrideDownloadOutcome Outcome,
        Exception Error)
    {
        /// <summary>The Aspera failure's kind when <see cref="Error"/> is a <see cref="PrideAsperaException"/>, else null.</summary>
        public PrideAsperaFailure? AsperaFailure => (Error as PrideAsperaException)?.Failure;
    }
}
