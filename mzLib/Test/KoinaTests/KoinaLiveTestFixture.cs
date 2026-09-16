using System;
using System.Net.Sockets;
using System.Threading.Tasks;
using NUnit.Framework;
using NUnit.Framework.Interfaces;
using NUnit.Framework.Internal;
using PredictionClients.Koina.Client;

namespace Test.KoinaTests
{
    /// <summary>
    /// Base class for test fixtures that call the live Koina inference server
    /// (https://koina.wilhelmlab.org). Two guards, because one is not enough:
    ///
    ///   * a [OneTimeSetUp] probe of the readiness endpoint, which skips the whole fixture when the
    ///     server is unreachable; and
    ///   * a [TearDown] that skips an individual test when the server was reachable and still could
    ///     not serve the request.
    ///
    /// The second exists because Koina's readiness endpoint answers for the server, not for the model,
    /// and a single model can be broken behind a healthy one. It has now been seen failing both ways
    /// a health check cannot report:
    ///
    ///   * mzLib #1241 -- Koina reported itself ready while one model's GPU was faulting, and every
    ///     call to that model came back 400 with a CUDA error in the body. It cost mzLib #1240 a red
    ///     check for the whole of its review.
    ///   * 2026-09-10 to 2026-09-11 -- /v2/models/AlphaPeptDeep_ccs_generic/ready answered 200 in
    ///     0.4 s while /v2/models/AlphaPeptDeep_ccs_generic/infer returned nothing at all, for
    ///     over a day. The session deadline expired and the four AlphaPeptDeep CCS tests failed
    ///     on every PR pushed in that window. Prosit_2023_intensity_XL_CMS2 and _CMS3 did the
    ///     same thing a day earlier, so it is not one model's problem.
    ///
    /// So the classification has to happen per call, and it has to cover a server that answers wrongly
    /// AND a server that does not answer.
    ///
    /// Deliberately narrower than <see cref="ExternalServiceTestHelper.RunAsync"/>, which skips on any
    /// HttpRequestException. A 400 meaning "your request is wrong" -- which is what a regression in how
    /// we build a request looks like -- still fails here.
    ///
    /// Derived fixtures still carry [Category("ExternalService")] and [Category("Koina")] so the
    /// CI category filters select them.
    /// </summary>
    public abstract class KoinaLiveTestFixture
    {
        private const string KoinaReadyUrl = "https://koina.wilhelmlab.org:443/v2/health/ready";

        /// <summary>
        /// Probes the live server once per fixture and skips the whole fixture if it is unreachable.
        /// </summary>
        /// <remarks>
        /// Virtual so that a fixture proving the guard works can suppress it. Hiding it with `new` does
        /// NOT: NUnit finds [OneTimeSetUp] by reflection and a `new` method leaves two distinct methods
        /// on the type, so BOTH run and the probe dials out anyway. An `override` leaves one.
        /// </remarks>
        [OneTimeSetUp]
        public virtual void EnsureKoinaReachable()
        {
            ExternalServiceTestHelper.EnsureReachable("Koina", KoinaReadyUrl);
        }

        /// <summary>
        /// Rewrites a failure caused by Koina being unable to serve the request into a skip.
        /// </summary>
        /// <remarks>
        /// A [TearDown] rather than a command-wrapper attribute because NUnit does not apply an
        /// IWrapTestMethod attribute declared on an abstract base class to the fixtures deriving from
        /// it -- tried, and the wrapper silently never fired, which is the failure mode that leaves
        /// every Koina test exactly as exposed as before with nothing to say so. [TearDown] inheritance
        /// is part of NUnit's documented contract, so this cannot go quietly wrong the same way.
        ///
        /// The Koina call happens deep inside model Predict(...) code, so there is no test body to
        /// wrap; the exception surfaces as the test's recorded result and is classified from there.
        /// </remarks>
        [TearDown]
        public void SkipWhenKoinaCouldNotServeTheRequest()
        {
            var result = TestExecutionContext.CurrentContext.CurrentResult;
            if (result.ResultState.Status != TestStatus.Failed) return;

            string? reason = ReasonKoinaCouldNotServeTheRequest(result);
            if (reason is null) return;

            string skipped =
                $"Skipping external-service test: {reason}. "
                + "This is a third-party availability problem, not a code failure.";
            TestContext.Progress.WriteLine(skipped);
            result.SetResult(ResultState.Ignored, skipped);
        }

        /// <summary>
        /// Why Koina could not serve this request, or null when the failure is ours to fix.
        /// </summary>
        /// <remarks>
        /// Two ways the server fails us, and they look nothing alike from here. It answered and could
        /// not run the model (a 400 carrying a CUDA fault), or it accepted the request and then never
        /// answered at all. The first is visible in the response body; the second has no body to read,
        /// which is why one classifier cannot cover both.
        /// </remarks>
        private static string? ReasonKoinaCouldNotServeTheRequest(ITestResult result)
        {
            string message = result.Message ?? string.Empty;

            // BOTH conditions, and neither alone is enough.
            //
            // The type name alone rewrites any failure whose message merely mentions the type -- most
            // obviously `Assert.Throws<KoinaServiceException>(...)` that threw nothing, which would be
            // reported Skipped when it is a genuine failure. Making the exception public in #1247 is
            // what put that test one keystroke away.
            //
            // The markers alone rewrite a test that quotes a CUDA error in an assertion message.
            //
            // Together they mean: the recorded failure names this exception AND carries a server fault
            // in its text. A Throws-nothing failure names the type but has no body, so it still fails.
            if (message.Contains(nameof(KoinaServiceException), StringComparison.Ordinal)
                && KoinaServiceException.IsServiceFault(message))
            {
                return $"Koina failed to run the model ({LineNamingTheFault(message)})";
            }

            // The same two-condition shape for the silent case, and for the same reason: a test that
            // asserted Throws<TaskCanceledException> and got nothing names the type in its message
            // too. What separates them is WHERE the abort was recorded -- inside the Koina client, or
            // in the test's own assertion.
            if (NamesATransportAbort(message) && AbortedInsideTheKoinaClient(result))
            {
                return $"Koina accepted the request and never answered it ({FirstLine(message)})";
            }

            return null;
        }

        /// <summary>
        /// Types that mean the call was cut off before any answer arrived.
        /// </summary>
        /// <remarks>
        /// Deliberately short, and it must stay that way, on the same admission test the response-body
        /// markers use: none of these is a plausible way for the server to report a request we built
        /// badly. A rejected request comes back as a 400 with a body, which is the branch above.
        ///
        /// IOException is NOT here even though the real failures nest one ("Unable to read data from
        /// the transport connection"). It is the one type on that chain a test could raise for an
        /// ordinary local reason -- reading a fixture file -- and the cancellation above it already
        /// names every case this exists for.
        /// </remarks>
        private static bool NamesATransportAbort(string message) =>
            message.Contains(nameof(TaskCanceledException), StringComparison.Ordinal)
            || message.Contains(nameof(OperationCanceledException), StringComparison.Ordinal)
            || message.Contains(nameof(TimeoutException), StringComparison.Ordinal)
            || message.Contains(nameof(SocketException), StringComparison.Ordinal);

        /// <summary>
        /// True when the recorded failure was raised inside the Koina client rather than by the test.
        /// </summary>
        /// <remarks>
        /// The stack trace, because it is the only thing that can tell a call that died from an
        /// assertion about one. An Assert.Throws&lt;TaskCanceledException&gt; that threw nothing is
        /// recorded against the test method and never enters this namespace; a call that actually
        /// aborted carries PredictionClients.Koina.Client.HTTP.InferenceRequest among its frames. The
        /// message is checked as well because Assert.DoesNotThrow reports the inner exception's trace
        /// inside the message rather than in StackTrace.
        ///
        /// What this deliberately does NOT separate: the session deadline expiring because Koina
        /// stalled, and it expiring because our own batch arithmetic underestimated the work. Both
        /// abort at the same line, so no amount of reading the failure tells them apart -- and
        /// skipping the second would be the mistake that got "out of memory" removed from
        /// KoinaServiceException.ServiceFaultMarkers. That guard is not weakened here, it is moved
        /// somewhere it works: KoinaModelBase.SessionDeadline is now one testable member, and
        /// KoinaModelDiscoveryTests.EveryConcreteModel_SessionDeadlineCoversItsOwnBenchmarkedWork
        /// pins it offline, so a batching regression fails a deterministic test in the required job
        /// rather than depending on a live one in the non-blocking job to notice.
        /// </remarks>
        private static bool AbortedInsideTheKoinaClient(ITestResult result) =>
            (result.StackTrace ?? string.Empty).Contains(KoinaClientNamespace, StringComparison.Ordinal)
            || (result.Message ?? string.Empty).Contains(KoinaClientNamespace, StringComparison.Ordinal);

        private const string KoinaClientNamespace = "PredictionClients.Koina";

        /// <summary>
        /// The first non-blank line of <paramref name="message"/>, which for a transport abort is the
        /// exception line itself.
        /// </summary>
        private static string FirstLine(string message)
        {
            foreach (string line in message.Split('\n'))
            {
                if (!string.IsNullOrWhiteSpace(line)) return Truncate(line.Trim());
            }
            return string.Empty;
        }

        /// <summary>
        /// The first line of <paramref name="message"/> that actually names the fault.
        /// </summary>
        /// <remarks>
        /// Taking the first line instead loses the diagnostic in exactly the case this exists for. An
        /// NUnit failure message is never valid JSON -- it is prefixed with the type name, or with
        /// "Expected: No Exception to be thrown" -- so ExtractServerError returns it unchanged, and its
        /// first line for an Assert.DoesNotThrow is the caller's user message. For the failure that
        /// motivated this PR that yielded "Charge 1 should be valid" and no mention of cuBLAS at all.
        /// </remarks>
        private static string LineNamingTheFault(string message)
        {
            foreach (string line in message.Split('\n'))
            {
                if (KoinaServiceException.IsServiceFault(line))
                {
                    return Truncate(line.Trim());
                }
            }
            return Truncate(message.Replace("\r", " ").Replace("\n", " ").Trim());
        }

        private static string Truncate(string text)
        {
            return text.Length <= 200 ? text : text[..200] + "...";
        }
    }
}
