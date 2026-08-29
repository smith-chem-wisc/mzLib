using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using NUnit.Framework.Interfaces;
using NUnit.Framework.Internal;
using PredictionClients.Koina.Client;

namespace Test.KoinaTests
{
    /// <summary>
    /// Covers the classification that mzLib #1241 was filed for: Koina answers 400 both for a request
    /// it will never accept and for a model it failed to run, so the status code cannot separate "not
    /// our bug, skip" from "our bug, fail". The body can.
    ///
    /// No network. The bodies here are the real ones, kept verbatim.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestKoinaServiceFaultClassification
    {
        /// <summary>
        /// The body that made #1240's CI red for its whole review, abbreviated only in the middle of
        /// the traceback.
        /// </summary>
        private const string RealCudaFault =
            "{\"error\":\"in ensemble 'Altimeter_2024_intensities', PyTorch execute failure: The "
            + "following operation failed in the TorchScript interpreter.\\nTraceback of TorchScript, "
            + "serialized code (most recent call last):\\n  File \\\"code/__torch__/models_spline.py\\\", "
            + "line 55, in forward\\n    out = torch.einsum(\\\"abc,bd->adc\\\", [inp, embed0])\\n"
            + "RuntimeError: CUDA error: CUBLAS_STATUS_EXECUTION_FAILED when calling "
            + "`cublasSgemmStridedBatched( handle, opa, opb, m, n, k, &alpha, a, lda, stridea, b, ldb, "
            + "strideb, &beta, c, ldc, stridec, num_batches)`\\n\"}";

        [Test]
        public void TheCudaFaultThatStartedThisIsClassifiedAsAServiceFault()
        {
            Assert.That(KoinaServiceException.IsServiceFault(RealCudaFault), Is.True);
            Assert.That(KoinaServiceException.ExtractServerError(RealCudaFault),
                Does.Contain("CUBLAS_STATUS_EXECUTION_FAILED"));
        }

        [TestCase("{\"error\":\"PyTorch execute failure: something\"}")]
        [TestCase("{\"error\":\"RuntimeError: CUDA error: device-side assert triggered\"}")]
        [TestCase("{\"error\":\"failed to load model 'Prosit_2019_intensity'\"}")]
        [TestCase("{\"error\":\"CUDA out of memory. Tried to allocate 2.00 GiB\"}")]
        public void ServerSideExecutionFailuresAreServiceFaults(string body)
        {
            Assert.That(KoinaServiceException.IsServiceFault(body), Is.True);
        }

        /// <summary>
        /// The half that matters more. A regression in how we build a request also arrives as a 400,
        /// and skipping on that would silence exactly the failure the suite exists to catch. Anything
        /// not recognised as a server-side execution failure has to fail.
        /// </summary>
        [TestCase("{\"error\":\"unexpected inference input 'peptide_sequenes' for model\"}")]
        [TestCase("{\"error\":\"expected 2 inputs but got 1\"}")]
        [TestCase("{\"error\":\"input 'peptide_sequences' batch size does not match other inputs\"}")]
        [TestCase("{\"error\":\"Model 'Prosit_2019_intensity' has no version '3'\"}")]
        [TestCase("")]
        [TestCase(null)]
        public void RequestRejectionsAreNotServiceFaults(string? body)
        {
            Assert.That(KoinaServiceException.IsServiceFault(body!), Is.False,
                "an unrecognised error must fail the test, not skip it");
        }

        /// <summary>
        /// A proxy or gateway can answer with HTML rather than Koina's JSON. That must be classified,
        /// not thrown on -- the body is the only evidence there is.
        /// </summary>
        [Test]
        public void ANonJsonBodyIsStillClassifiedRatherThanThrowing()
        {
            const string html = "<html><body><h1>502 Bad Gateway</h1></body></html>";

            Assert.That(() => KoinaServiceException.IsServiceFault(html), Throws.Nothing);
            Assert.That(KoinaServiceException.IsServiceFault(html), Is.False);
            Assert.That(KoinaServiceException.ExtractServerError(html), Is.EqualTo(html),
                "with no error field, the raw body is the best description available");
        }

        [Test]
        public void ExistingCatchClausesStillSeeItAsAnHttpRequestException()
        {
            var e = new KoinaServiceException("Request failed with status 400", "CUDA error");

            Assert.That(e, Is.InstanceOf<System.Net.Http.HttpRequestException>(),
                "callers already catching HttpRequestException must keep working unchanged");
            Assert.That(e.ServerError, Is.EqualTo("CUDA error"));
        }
    }

    /// <summary>
    /// Proves the inherited guard actually fires. This is the part most likely to be silently wrong,
    /// and was: the first version used an IWrapTestMethod attribute on the base class, NUnit did not
    /// apply it to derived fixtures, and the wrapper never ran -- which would have left every Koina
    /// test exactly as exposed as before with nothing to say so. This test is what caught that.
    ///
    /// Derives from <see cref="KoinaLiveTestFixture"/> the same way the real fixtures do, but throws
    /// the fault itself rather than calling Koina, so it needs no network and runs in the ordinary
    /// suite.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestKoinaServiceFaultIsSkippedNotFailed : KoinaLiveTestFixture
    {
        /// <summary>
        /// The base class's [OneTimeSetUp] probes the live server, which this fixture must not do.
        /// </summary>
        [OneTimeSetUp]
        public new void EnsureKoinaReachable() { }

        [Test]
        public void AServiceFaultSkipsTheTest()
        {
            throw new KoinaServiceException(
                "Request failed with status 400 Bad Request",
                "PyTorch execute failure: RuntimeError: CUDA error: CUBLAS_STATUS_EXECUTION_FAILED");
        }

        /// <summary>
        /// The wrapper must not swallow anything else. If this one ever reports Skipped, the wrapper
        /// has become the over-skipping the issue warned against.
        /// </summary>
        [Test]
        public void AnOrdinaryFailureStillFails()
        {
            Assert.That(1, Is.EqualTo(1));
        }

        /// <summary>
        /// Observes the finished results rather than doing it in a [TearDown]: NUnit runs TearDowns
        /// derived-first, base-last, so a TearDown in this fixture would look at the result BEFORE the
        /// base class's guard had rewritten it. That is what the first version of this test did, and it
        /// reported a failure that was its own.
        /// </summary>
        [OneTimeTearDown]
        public void TheGuardMustHaveSkippedTheFaultAndLeftTheOtherAlone()
        {
            var children = TestExecutionContext.CurrentContext.CurrentResult.Children.ToList();

            var fault = children.Single(c => c.Name == nameof(AServiceFaultSkipsTheTest));
            Assert.That(fault.ResultState.Status, Is.EqualTo(TestStatus.Skipped),
                "a Koina service fault must be skipped, not failed -- the guard did not fire");
            Assert.That(fault.Message, Does.Contain("third-party availability problem"));

            var ordinary = children.Single(c => c.Name == nameof(AnOrdinaryFailureStillFails));
            Assert.That(ordinary.ResultState.Status, Is.EqualTo(TestStatus.Passed),
                "the guard must not touch anything that is not a Koina service fault");
        }
    }
}
