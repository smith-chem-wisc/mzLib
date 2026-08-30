using System;
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
    ///   * a [TearDown] that skips an individual test when the server answered but failed to RUN the
    ///     model.
    ///
    /// The second exists because of mzLib #1241. Koina reported itself ready while one model's GPU was
    /// faulting, and every call to that model came back 400 with a CUDA error in the body. A health
    /// endpoint cannot report that, so a fixture-level probe cannot catch it -- the classification has
    /// to happen per call. It cost mzLib #1240 a red check for the whole of its review.
    ///
    /// Deliberately narrower than <see cref="ExternalServiceTestHelper.RunAsync"/>, which skips on any
    /// HttpRequestException. Only a <see cref="KoinaServiceException"/> is skipped here, so a 400
    /// meaning "your request is wrong" -- which is what a regression in how we build a request looks
    /// like -- still fails.
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
        /// Rewrites a failure caused by Koina failing to run the model into a skip.
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
        public void SkipWhenKoinaFailedToRunTheModel()
        {
            var result = TestExecutionContext.CurrentContext.CurrentResult;
            if (result.ResultState.Status != TestStatus.Failed) return;

            // BOTH conditions, and neither alone is enough.
            //
            // The type name alone rewrites any failure whose message merely mentions the type -- most
            // obviously `Assert.Throws<KoinaServiceException>(...)` that threw nothing, which would be
            // reported Skipped when it is a genuine failure. Making the exception public in this PR is
            // what puts that test one keystroke away.
            //
            // The markers alone rewrite a test that quotes a CUDA error in an assertion message.
            //
            // Together they mean: the recorded failure names this exception AND carries a server fault
            // in its text. A Throws-nothing failure names the type but has no body, so it still fails.
            string message = result.Message ?? string.Empty;
            if (!message.Contains(nameof(KoinaServiceException), StringComparison.Ordinal)) return;
            if (!KoinaServiceException.IsServiceFault(message)) return;

            string skipped =
                "Skipping external-service test: Koina failed to run the model "
                + $"({LineNamingTheFault(message)}). "
                + "This is a third-party availability problem, not a code failure.";
            TestContext.Progress.WriteLine(skipped);
            result.SetResult(ResultState.Ignored, skipped);
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
