using System;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// How <see cref="PrideArchiveClient.DownloadFileViaAsperaAsync"/> runs IBM Aspera's <c>ascp</c> to fetch a PRIDE
    /// file from EBI's Aspera server.
    /// </summary>
    /// <remarks>
    /// mzLib does not ship ascp. EBI (KB0011565) requires ascp 4.4.4 or newer with encryption on: IBM Aspera for
    /// Desktop installs 4.4.8 at <c>C:\Program Files\IBM Aspera\transferd\bin\ascp.exe</c>, and on Linux IBM's
    /// Transfer SDK installs it unattended (<c>gem install aspera-cli; ascli config ascp install</c>). Aspera Connect's
    /// 4.4.3 failed to authenticate against EBI, and so did the 3.7.3 that pridepy redistributes for Linux.
    /// </remarks>
    public sealed class PrideAsperaOptions
    {
        /// <summary>The full path of the <c>ascp</c> executable. Required; mzLib does not search for one.</summary>
        public string AscpPath { get; init; }

        /// <summary>
        /// The private key ascp logs in with. Null (the default) uses EBI's public key for anonymous PRIDE downloads,
        /// which mzLib bundles; set it only to use another copy of that key.
        /// </summary>
        public string KeyPath { get; init; }

        /// <summary>ascp's target rate, its <c>-l</c> argument: digits, optionally followed by k, m or g. "300m" by default.</summary>
        public string MaxRate { get; init; } = "300m";

        /// <summary>
        /// How long the file may receive no bytes before ascp is stopped and the attempt fails as
        /// <see cref="PrideAsperaFailure.Stalled"/>. An inactivity deadline, not a duration one: a slow transfer that
        /// keeps receiving is never stopped. Five minutes by default.
        /// </summary>
        /// <remarks>
        /// ascp can hang without exiting: on 2026-10-29 one session sat 9.3 hours with no bytes while EBI refused SSH
        /// connections. Five minutes is the value PXReprise's batch has run since then; ascp's own login takes
        /// seconds, so a session that has received nothing for five minutes is not about to.
        /// </remarks>
        public TimeSpan StallTimeout { get; init; } = TimeSpan.FromMinutes(5);
    }
}
