using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Reflection;
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
        internal const string RealCudaFault =
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

        /// <summary>
        /// The decision itself, which is the whole point of the type: a service fault becomes a
        /// KoinaServiceException, and anything else stays the plain exception it always was. Tested
        /// through the factory rather than through InferenceRequest, because provoking a live server
        /// into failing on demand is not something a test can do.
        /// </summary>
        [Test]
        public void AServiceFaultBecomesAKoinaServiceExceptionCarryingTheServerError()
        {
            var thrown = KoinaServiceException.ForFailedResponse(
                400, "Bad Request", TestKoinaServiceFaultClassification.RealCudaFault);

            Assert.That(thrown, Is.InstanceOf<KoinaServiceException>());
            Assert.That(((KoinaServiceException)thrown).ServerError,
                Does.Contain("CUBLAS_STATUS_EXECUTION_FAILED"));
            Assert.That(thrown.Message, Does.Contain("400").And.Contain("Bad Request"),
                "the status still has to reach the reader");
        }

        [Test]
        public void ARequestRejectionStaysAPlainHttpRequestException()
        {
            var thrown = KoinaServiceException.ForFailedResponse(
                400, "Bad Request", "{\"error\":\"expected 2 inputs but got 1\"}");

            Assert.That(thrown, Is.Not.InstanceOf<KoinaServiceException>(),
                "our own malformed request must not be dressed up as the server's fault");
            Assert.That(thrown, Is.InstanceOf<System.Net.Http.HttpRequestException>());
        }

        /// <summary>
        /// A 500 with no body at all still has to produce something a caller can read. It is not a
        /// service fault by this classifier -- nothing in the body says the model failed, because there
        /// is no body -- and that is the safe direction.
        /// </summary>
        [Test]
        public void AnEmptyBodyStillProducesAReadableFailure()
        {
            var thrown = KoinaServiceException.ForFailedResponse(500, "Internal Server Error", "");

            Assert.That(thrown, Is.Not.InstanceOf<KoinaServiceException>());
            Assert.That(thrown.Message, Does.Contain("500").And.Contain("Internal Server Error"));
        }

        /// <summary>
        /// JSON that is well formed but not the shape Koina documents. Both of these parse, so the
        /// JsonReaderException path does not catch them; the body itself is returned instead of a cast
        /// blowing up on a token that is not a string.
        /// </summary>
        [TestCase("[{\"error\":\"CUDA error\"}]", Description = "a JSON array, not an object")]
        [TestCase("{\"error\":{\"code\":500}}", Description = "an error field that is not a string")]
        [TestCase("{\"detail\":\"CUDA error\"}", Description = "no error field at all")]
        public void WellFormedJsonOfTheWrongShapeFallsBackToTheRawBody(string body)
        {
            Assert.That(() => KoinaServiceException.ExtractServerError(body), Throws.Nothing);
            Assert.That(KoinaServiceException.ExtractServerError(body), Is.EqualTo(body));
        }

        /// <summary>
        /// A GPU out-of-memory is NOT a service fault, because it is exactly what an oversized request
        /// looks like from the server side -- too large a batch, too many peptides in one call, a
        /// pathological sequence length. The body says out of memory and the fault is ours.
        ///
        /// That also makes it the plausible regression here: batch size gets tuned, Koina OOMs, and CI
        /// would go green with a skip. Co-occurrence with a cuda marker does not separate the two,
        /// since a request-induced OOM reports itself the same way.
        /// </summary>
        [TestCase("{\"error\":\"CUDA out of memory. Tried to allocate 2.00 GiB\"}")]
        [TestCase("{\"error\":\"RuntimeError: CUDA out of memory\"}")]
        public void AnOutOfMemoryIsOurProblemNotTheServersAndMustNotSkip(string body)
        {
            Assert.That(KoinaServiceException.IsServiceFault(body), Is.False,
                "an oversized request OOMs the GPU too, so this has to fail rather than skip");
        }

        /// <summary>
        /// Triton names the model in the middle -- "model 'X' is not ready" -- so a marker requiring
        /// those four words adjacent would never fire.
        /// </summary>
        [Test]
        public void AnUnreadyModelIsRecognisedWithTheModelNameInTheMiddle()
        {
            Assert.That(
                KoinaServiceException.IsServiceFault(
                    "{\"error\":\"model 'Altimeter_2024_intensities' is not ready\"}"),
                Is.True);
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
        /// Suppresses the base class's live probe, which this fixture must not perform.
        /// </summary>
        /// <remarks>
        /// `override`, not `new`. Hiding it leaves two [OneTimeSetUp] methods on the type and NUnit runs
        /// BOTH, so the probe dialled out anyway -- and this fixture carries no [Category], so
        /// `cat != ExternalService` selects it into the required job that #1140 exists to keep off the
        /// network. Worse, EnsureReachable calls Assert.Ignore, so with Koina down the fixture proving
        /// the guard works was itself skipped: the guard unproven, with nothing saying so.
        /// </remarks>
        [OneTimeSetUp]
        public override void EnsureKoinaReachable() { }

        [Test]
        public void AServiceFaultSkipsTheTest()
        {
            // Built through the factory rather than by hand, because the message shape matters: the
            // guard now requires the recorded failure to carry a server fault in its TEXT, and it is
            // ForFailedResponse that puts the response body there. Constructing the exception directly
            // with a body-free message tested a shape production never produces.
            throw KoinaServiceException.ForFailedResponse(
                400, "Bad Request", TestKoinaServiceFaultClassification.RealCudaFault);
        }

        /// <summary>
        /// The shape that motivated the PR. TestAltimeterChargeStateBoundaries wraps the live call in
        /// Assert.DoesNotThrow with a user message, so the recorded failure begins "Expected: No
        /// Exception to be thrown" and its FIRST line is the user message, not the server's.
        ///
        /// Reporting the first line lost the diagnostic exactly here: the skip read "(  Charge 1 should
        /// be valid)" and cuBLAS appeared nowhere. The guard now reports the first line that names the
        /// fault, so the log says what Koina actually did.
        /// </summary>
        [Test]
        public void ADoesNotThrowFailureStillReportsTheServersError()
        {
            Assert.DoesNotThrow(() => throw KoinaServiceException.ForFailedResponse(
                400, "Bad Request", TestKoinaServiceFaultClassification.RealCudaFault),
                "Charge 1 should be valid");
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
        /// A test that expected this exception and did not get one is a GENUINE failure and must stay
        /// one. Its NUnit message names the type -- "Expected: &lt;...KoinaServiceException&gt;" -- so a
        /// guard keyed on the type name alone would report it Skipped.
        ///
        /// Latent rather than live when this was written: every Throws in the live fixtures is
        /// Assert.Throws&lt;ArgumentException&gt;. But this PR makes KoinaServiceException public, and
        /// Assert.Throws&lt;KoinaServiceException&gt; against InferenceRequest is the next test someone
        /// writes.
        /// </summary>
        [Test]
        public void AFailedThrowsAssertionNamingTheExceptionIsStillAFailure()
        {
            var recorded = Assert.Throws<AssertionException>(() =>
                Assert.That(() => { }, Throws.InstanceOf<KoinaServiceException>()));

            Assert.That(recorded!.Message, Does.Contain(nameof(KoinaServiceException)),
                "premise: the message names the type, which is what could mislead the guard");
            Assert.That(KoinaServiceException.IsServiceFault(recorded.Message), Is.False,
                "and it carries no server fault, which is what keeps the guard off it");
        }

        /// <summary>
        /// Hiding the base [OneTimeSetUp] with `new` leaves TWO of them on the type and NUnit runs both,
        /// so the live probe fired from a fixture documented as needing no network. Reflection is how
        /// that is visible, and it is the same lookup NUnit does.
        /// </summary>
        [Test]
        public void TheLiveProbeIsSuppressedRatherThanHidden()
        {
            // Flags rather than the parameterless GetMethods(), which sees public methods only: were the
            // base probe ever narrowed to `protected virtual`, that overload would return zero and the
            // count assertion below would fail with a message about `new` versus `override` describing
            // something that did not happen.
            var probes = typeof(TestKoinaServiceFaultIsSkippedNotFailed)
                .GetMethods(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.Static)
                .Where(m => m.Name == nameof(KoinaLiveTestFixture.EnsureKoinaReachable))
                .ToList();

            Assert.That(probes, Has.Count.EqualTo(1),
                "two of them means `new` rather than `override`, and NUnit would run the base probe too");
            Assert.That(probes[0].DeclaringType, Is.EqualTo(typeof(TestKoinaServiceFaultIsSkippedNotFailed)));
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

            // SingleOrDefault, not Single: a --filter that selects one test in this fixture still runs
            // this [OneTimeTearDown], and a Single over an absent sibling would throw from the teardown.
            // The person debugging one of these would then be shown a failure that is the observer's
            // rather than theirs -- the same trap the doc comment above records this test falling into.
            // Each assertion therefore runs only over the children the run actually produced.
            ITestResult Child(string name) => children.SingleOrDefault(c => c.Name == name);

            var fault = Child(nameof(AServiceFaultSkipsTheTest));
            if (fault is not null)
            {
                Assert.That(fault.ResultState.Status, Is.EqualTo(TestStatus.Skipped),
                    "a Koina service fault must be skipped, not failed -- the guard did not fire");
                Assert.That(fault.Message, Does.Contain("third-party availability problem"));
            }

            var ordinary = Child(nameof(AnOrdinaryFailureStillFails));
            if (ordinary is not null)
            {
                Assert.That(ordinary.ResultState.Status, Is.EqualTo(TestStatus.Passed),
                    "the guard must not touch anything that is not a Koina service fault");
            }

            var wrapped = Child(nameof(ADoesNotThrowFailureStillReportsTheServersError));
            if (wrapped is not null)
            {
                Assert.That(wrapped.ResultState.Status, Is.EqualTo(TestStatus.Skipped));
                Assert.That(wrapped.Message, Does.Contain("PyTorch execute failure"),
                    "the server's own error has to reach the log, not the caller's user message");
                Assert.That(wrapped.Message, Does.Not.Contain("Charge 1 should be valid"));
            }

            var throwsNothing = Child(nameof(AFailedThrowsAssertionNamingTheExceptionIsStillAFailure));
            if (throwsNothing is not null)
            {
                Assert.That(throwsNothing.ResultState.Status, Is.EqualTo(TestStatus.Passed));
            }
        }
    }
}
