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
    /// <param name="Elapsed">From the start of the transfer to its end, not counting waits for an Aspera session.</param>
    /// <param name="Outcome">How it ended.</param>
    /// <param name="Error">Why it did not succeed; null when <see cref="PrideDownloadOutcome.Accepted"/>.</param>
    public sealed record PrideDownloadTry(string FileName, int Try, PrideDownloadRoute Route, string Host,
        long FileSizeBytes, TimeSpan Elapsed, PrideDownloadOutcome Outcome, Exception Error)
    {
        /// <summary>The Aspera failure's kind when <see cref="Error"/> is a <see cref="PrideAsperaException"/>, else null.</summary>
        public PrideAsperaFailure? AsperaFailure => (Error as PrideAsperaException)?.Failure;
    }
}
