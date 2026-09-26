using System.Collections.Generic;
using System.Threading;
using NUnit.Framework;
using PredictionClients.Koina.Client;

namespace Test.KoinaTests
{
    /// <summary>
    /// No-network tests for the Koina HTTP wrapper. The wrapper now exposes a process-wide shared
    /// client with per-request cancellation; the request methods themselves are exercised by the
    /// Koina-category integration tests.
    /// </summary>
    [TestFixture]
    public class HttpClientUnitTests
    {
        [Test]
        public void ModelsUrl_PointsAtKoina()
        {
            Assert.That(HTTP.ModelsURL, Does.StartWith("https://koina.wilhelmlab.org"));
            Assert.That(HTTP.ModelsURL, Does.EndWith("/models/"));
        }

        [Test]
        public void InferenceRequest_AlreadyCanceledToken_Throws()
        {
            // A pre-canceled token short-circuits before any network I/O, so this stays offline.
            using var cts = new CancellationTokenSource();
            cts.Cancel();

            Assert.CatchAsync<System.OperationCanceledException>(
                async () => await HTTP.InferenceRequest("any_model", new Dictionary<string, object>(), cts.Token));
        }
    }
}
