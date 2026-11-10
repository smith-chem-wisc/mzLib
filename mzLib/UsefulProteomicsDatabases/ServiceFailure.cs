using System;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Tells a repository client's failure that is worth trying again later apart from one that will not change,
    /// for the exceptions <see cref="PrideArchiveClient"/> throws.
    /// </summary>
    public static class ServiceFailure
    {
        /// <summary>The <see cref="Exception.Data"/> key the retry loop marks a transient failure it gave up on with.</summary>
        internal const string TransientKey = "ServiceFailure.Transient";

        /// <summary>
        /// Whether <paramref name="exception"/> is a failure the client judged transient and stopped retrying only
        /// because its retries ran out: the service timing out, a body that stalled, dropped or came up short, a
        /// 408, a 429, a 5xx, or PRIDE's 403 on its FTP host. Asking again later may well succeed.
        /// </summary>
        /// <remarks>
        /// Read from a mark the client's retry loop leaves on the exception it gives up on, not from the exception's
        /// shape, because the shape cannot tell:
        /// <list type="bullet">
        /// <item>a 403 is transient on PRIDE's FTP host (EBI's rate limit) and a refusal everywhere else, and the
        /// exception does not carry the host it came from;</item>
        /// <item>an <see cref="System.Net.Http.HttpRequestException"/> with no status is also what the client throws
        /// for a few answers it never retries (a PROXI response with no spectrum, more pages than
        /// <see cref="PrideArchiveClient.MaxPages"/>).</item>
        /// </list>
        /// False for everything else: an <c>MzLibException</c> (the service answered and broke its contract, or a
        /// download failed its checks on every attempt, which is a data error, not an outage), an
        /// <see cref="ArgumentException"/>, a refusal such as 404, the caller's cancellation, a local disk error,
        /// an exception that did not come from a mzLib client, and null.
        /// <para>
        /// Only <paramref name="exception"/> itself is read, never its <see cref="Exception.InnerException"/>: the
        /// <c>MzLibException</c> for a damaged download wraps its last attempt, which was transient, but the
        /// verdict on the download is that it is wrong.
        /// </para>
        /// </remarks>
        /// <param name="exception">An exception thrown by a mzLib repository client; null gives false.</param>
        public static bool IsTransient(Exception exception) => exception?.Data[TransientKey] is true;

        /// <summary>Marks <paramref name="exception"/> as a transient failure the retry loop gave up on.</summary>
        internal static void MarkTransient(Exception exception) => exception.Data[TransientKey] = true;
    }
}
