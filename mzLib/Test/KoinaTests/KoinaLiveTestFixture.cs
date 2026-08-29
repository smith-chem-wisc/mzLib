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

        [OneTimeSetUp]
        public void EnsureKoinaReachable()
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

            // The recorded message carries both the exception type name and the server's response
            // body. Either alone identifies the fault; requiring the type name keeps this from firing
            // on a test that merely mentions CUDA in an assertion message.
            string message = result.Message ?? string.Empty;
            if (!message.Contains(nameof(KoinaServiceException), StringComparison.Ordinal)) return;

            string skipped =
                "Skipping external-service test: Koina failed to run the model "
                + $"({FirstLine(KoinaServiceException.ExtractServerError(message))}). "
                + "This is a third-party availability problem, not a code failure.";
            TestContext.Progress.WriteLine(skipped);
            result.SetResult(ResultState.Ignored, skipped);
        }

        private static string FirstLine(string text)
        {
            int newline = text.IndexOfAny(['\r', '\n']);
            string line = newline < 0 ? text : text[..newline];
            return line.Length <= 200 ? line : line[..200] + "...";
        }
    }
}
