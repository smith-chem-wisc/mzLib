using NUnit.Framework;

namespace Test.KoinaTests
{
    /// <summary>
    /// Base class for test fixtures that call the live Koina inference server
    /// (https://koina.wilhelmlab.org). A single [OneTimeSetUp] probe checks the server's
    /// readiness endpoint once per fixture; if Koina is unreachable the whole fixture is
    /// Skipped (with a reason on the CI log) rather than failing - see
    /// <see cref="ExternalServiceTestHelper"/>. Koina calls are made deep inside model
    /// Predict(...) code, so a fixture-level probe is the least invasive way to distinguish
    /// "Koina is down" (skip) from a genuine regression (fail).
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
    }
}
