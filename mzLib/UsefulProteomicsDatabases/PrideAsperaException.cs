using System;
using System.Net.Http;

namespace UsefulProteomicsDatabases
{
    /// <summary>What went wrong with one Aspera download, as <see cref="PrideAsperaException.Failure"/> reports it.</summary>
    /// <remarks>
    /// The kinds a caller acts on differently. <see cref="LoginRefused"/> and <see cref="CannotConnect"/> are failures of
    /// the route, not of the file: every file fails the same way until EBI recovers, which took about three days for the
    /// refusal of 2026-10-30. The others belong to one session, and another session may well succeed.
    /// </remarks>
    public enum PrideAsperaFailure
    {
        /// <summary>
        /// EBI refused the key: ascp's <c>Session Stop (Error: failed to authenticate)</c>. Seen 2,272 times in
        /// PXReprise's batch from 2026-10-30 to 2026-11-03, against a key that never changed.
        /// </summary>
        LoginRefused,

        /// <summary>
        /// ascp could not reach the server: <c>Unable to connect via SSH</c>, <c>Failed to open TCP connection for
        /// SSH</c>, <c>Client unable to connect to server (check UDP port and firewall)</c>, or <c>Session initiation
        /// failed</c>.
        /// </summary>
        CannotConnect,

        /// <summary>No bytes arrived for <see cref="PrideAsperaOptions.StallTimeout"/>, so ascp was stopped.</summary>
        Stalled,

        /// <summary>
        /// ascp finished, but the file failed the check against PRIDE's checksum row: the wrong size, the wrong MD5, or
        /// EBI's error text written into it. The copy was discarded.
        /// </summary>
        DamagedCopy,

        /// <summary>
        /// ascp failed in any other way, for example <c>Session data transfer timeout</c>, or exited without writing
        /// the file.
        /// </summary>
        AscpFailed
    }

    /// <summary>
    /// One Aspera download of a PRIDE file that failed in a way worth trying again later; see <see cref="Failure"/> for
    /// which. <see cref="ServiceFailure.IsTransient"/> is true for every one.
    /// </summary>
    /// <remarks>
    /// An <see cref="HttpRequestException"/> with no status, the type mzLib's clients use for a transport failure, so a
    /// caller and <c>ExternalServiceTestHelper</c> treat it as they treat an HTTPS outage. The message names the file and
    /// <see cref="PrideArchiveExtensions.PrideAsperaHost"/>, never the source path. The
    /// <see cref="Exception.InnerException"/> carries what ascp said, or the check that failed.
    /// </remarks>
    public sealed class PrideAsperaException : HttpRequestException
    {
        internal PrideAsperaException(string message, PrideAsperaFailure failure, int? ascpExitCode, Exception innerException)
            : base(message, innerException)
        {
            Failure = failure;
            AscpExitCode = ascpExitCode;
            ServiceFailure.MarkTransient(this);
        }

        /// <summary>What went wrong.</summary>
        public PrideAsperaFailure Failure { get; }

        /// <summary>
        /// ascp's exit code when it exited on its own (0 for a <see cref="PrideAsperaFailure.DamagedCopy"/>), else null:
        /// it was stopped.
        /// </summary>
        public int? AscpExitCode { get; }
    }
}
