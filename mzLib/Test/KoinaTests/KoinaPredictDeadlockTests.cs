using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CCSModels;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;

namespace Test.KoinaTests
{
    /// <summary>
    /// Regression guard for the synchronous Predict() wrappers across every model family. Each must
    /// offload its async work via Task.Run so it cannot deadlock when called from a single-threaded
    /// SynchronizationContext (WinForms / WPF / ASP.NET non-Core).
    ///
    /// The fake transports yield asynchronously (await Task.Yield()) so that a naive
    /// .GetAwaiter().GetResult() would post a continuation onto the captured, blocked context and hang.
    /// Correct code runs the continuations on the ThreadPool (no context) and returns. If the Task.Run
    /// offload is ever removed from any Predict(), that family's test fails via the Join timeout.
    /// </summary>
    [TestFixture]
    public class KoinaPredictDeadlockTests
    {
        private const string IntensityJson =
            "{\"outputs\":[" +
            "{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[2],\"data\":[\"b1+1\",\"y1+1\"]}," +
            "{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[2],\"data\":[100.0,200.0]}," +
            "{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[2],\"data\":[0.5,0.6]}]}";

        [Test]
        public void FragmentIntensity_Predict_UnderSingleThreadedContext_DoesNotDeadlock() =>
            AssertCompletesUnderSingleThreadedContext(() => new YieldingFragmentModel()
                .Predict(new List<FragmentIntensityPredictionInput> { new("PEPTIDEK", 2, 30, null, null) }));

        [Test]
        public void RetentionTime_Predict_UnderSingleThreadedContext_DoesNotDeadlock() =>
            AssertCompletesUnderSingleThreadedContext(() => new YieldingRtModel()
                .Predict(new List<RetentionTimePredictionInput> { new("PEPTIDEK") }));

        [Test]
        public void Ccs_Predict_UnderSingleThreadedContext_DoesNotDeadlock() =>
            AssertCompletesUnderSingleThreadedContext(() => new YieldingCcsModel()
                .Predict(new List<CCSPredictionInput> { new("PEPTIDEK", 2) }));

        [Test]
        public void Crosslink_Predict_UnderSingleThreadedContext_DoesNotDeadlock() =>
            AssertCompletesUnderSingleThreadedContext(() => new YieldingCrosslinkModel()
                .Predict(new List<CrosslinkIntensityPredictionInput> { new("PEPTIDEK[UNIMOD:1896]", "ACDEK[UNIMOD:1896]", 2, 30) }));

        [Test]
        public void Detectability_Predict_UnderSingleThreadedContext_DoesNotDeadlock() =>
            AssertCompletesUnderSingleThreadedContext(() => new YieldingDetectabilityModel()
                .Predict(new List<DetectabilityPredictionInput> { new("PEPTIDEK") }));

        // Runs predict() on a dedicated thread whose SynchronizationContext captures posted continuations
        // but never pumps them (the thread is blocked inside Predict). Correct code offloads to the
        // ThreadPool and returns; a direct GetResult() waits forever for a continuation needing this thread.
        private static void AssertCompletesUnderSingleThreadedContext(Action predict)
        {
            Exception? captured = null;
            var thread = new Thread(() =>
            {
                SynchronizationContext.SetSynchronizationContext(new QueueingSynchronizationContext());
                try { predict(); }
                catch (Exception ex) { captured = ex; }
            }) { IsBackground = true };

            thread.Start();
            bool completed = thread.Join(TimeSpan.FromSeconds(15));

            Assert.That(completed, Is.True,
                "Predict did not return under a single-threaded SynchronizationContext — likely a .GetAwaiter().GetResult() deadlock.");
            if (captured != null) throw captured;
        }

        // Captures posted continuations on a queue that is never drained (the owning thread blocks in
        // Predict). Reproduces the WinForms/WPF condition where the UI thread cannot pump while blocked.
        private sealed class QueueingSynchronizationContext : SynchronizationContext
        {
            private readonly BlockingCollection<(SendOrPostCallback d, object? s)> _queue = new();
            public override void Post(SendOrPostCallback d, object? state) => _queue.Add((d, state));
            public override void Send(SendOrPostCallback d, object? state) => d(state);
        }

        // Fake transports yield before returning a canned response so the awaited continuation is posted
        // to the captured context (rather than completing synchronously, which would not exercise the bug).
        private sealed class YieldingFragmentModel : Prosit2019Intensity
        {
            protected override async Task<string> SendInferenceRequestAsync(string m, Dictionary<string, object> r, CancellationToken ct)
            { await Task.Yield(); return IntensityJson; }
        }

        private sealed class YieldingRtModel : Prosit2019iRT
        {
            protected override async Task<string> SendInferenceRequestAsync(string m, Dictionary<string, object> r, CancellationToken ct)
            { await Task.Yield(); return "{\"outputs\":[{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[1],\"data\":[42.5]}]}"; }
        }

        private sealed class YieldingCcsModel : IM2Deep
        {
            protected override async Task<string> SendInferenceRequestAsync(string m, Dictionary<string, object> r, CancellationToken ct)
            { await Task.Yield(); return "{\"outputs\":[{\"name\":\"ccs\",\"datatype\":\"FP32\",\"shape\":[1],\"data\":[321.5]}]}"; }
        }

        private sealed class YieldingCrosslinkModel : Prosit2023IntensityXLCMS2
        {
            protected override async Task<string> SendInferenceRequestAsync(string m, Dictionary<string, object> r, CancellationToken ct)
            { await Task.Yield(); return IntensityJson; }
        }

        private sealed class YieldingDetectabilityModel : PFly2024FineTuned
        {
            protected override async Task<string> SendInferenceRequestAsync(string m, Dictionary<string, object> r, CancellationToken ct)
            { await Task.Yield(); return "{\"outputs\":[{\"name\":\"detectability\",\"datatype\":\"FP32\",\"shape\":[4],\"data\":[0.1,0.2,0.3,0.4]}]}"; }
        }
    }
}
