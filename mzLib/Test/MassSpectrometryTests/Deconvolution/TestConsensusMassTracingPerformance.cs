using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MzLibUtil;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Step-1 PERFORMANCE / SCALING driver for the consensus mass-tracing engines
    /// (plan: E:\CodeReview\PR1069_test_and_perf_plan.md). [Explicit] so it never
    /// runs in CI timing: it allocates large synthetic inputs and prints stage
    /// timings plus a doubling-scaling ratio.
    ///
    /// Purpose: (a) get real isolated numbers for BuildTraces / Correct /
    /// BuildFeatures (the Velos E2E numbers fold in Classic decon, which dominates);
    /// (b) catch an accidental regression toward O(n^2) in BuildTraces' open-trace
    /// inner loop.
    ///
    /// Run: dotnet test --filter FullyQualifiedName~TestConsensusMassTracingPerformance
    /// </summary>
    [TestFixture]
    [Explicit("Performance/scaling driver; allocates large synthetic data and prints timings.")]
    [ExcludeFromCodeCoverage]
    public class TestConsensusMassTracingPerformance
    {
        private const double PassBToleranceDa = 1.5;
        private const int MaxGap = 1;
        private const double CrossChargeMassPpm = 10.0;

        [Test]
        public void StageTimings_RepresentativeTopDownSize()
        {
            RunOnce(nScans: 2000, activeSpeciesPerScan: 60, chargesPerSpecies: 3, label: "TD ~2000 scans");
        }

        [Test]
        public void Scaling_BuildTraces_IsSubQuadratic()
        {
            // Open-trace count is held constant (density fixed), so doubling the scan
            // count should grow BuildTraces ~linearly. A quadratic regression in the
            // open-trace inner loop would push the ratio toward 4.
            var small = RunOnce(nScans: 1000, activeSpeciesPerScan: 60, chargesPerSpecies: 3, label: "N=1000");
            var large = RunOnce(nScans: 2000, activeSpeciesPerScan: 60, chargesPerSpecies: 3, label: "N=2000");

            double ratio = large.BuildTracesMs / Math.Max(1.0, small.BuildTracesMs);
            TestContext.Out.WriteLine($"BuildTraces scaling  T(2N)/T(N) = {ratio:F2}   (linear≈2, quadratic≈4)");

            Assert.That(ratio, Is.LessThan(3.0),
                "BuildTraces scaled worse than ~linear when scans doubled -- possible O(n^2) regression in the open-trace loop");
        }

        private sealed record StageResult(
            int Envelopes, int Traces, int Features, int MultiCharge,
            double BuildTracesMs, double CorrectMs, double BuildFeaturesMs);

        private static StageResult RunOnce(int nScans, int activeSpeciesPerScan, int chargesPerSpecies, string label)
        {
            var (scans, perScan, totalEnv) = Generate(nScans, activeSpeciesPerScan, chargesPerSpecies);

            var sw = Stopwatch.StartNew();
            var traces = MassTraceBuilder.BuildTraces(scans, perScan, PassBToleranceDa, MaxGap);
            sw.Stop(); double buildMs = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();
            sw.Stop(); double corrMs = sw.Elapsed.TotalMilliseconds;

            sw.Restart();
            var features = MassFeatureBuilder.BuildFeatures(corrected, CrossChargeMassPpm);
            sw.Stop(); double featMs = sw.Elapsed.TotalMilliseconds;

            int multi = features.Count(f => f.ChargeCount >= 2);
            TestContext.Out.WriteLine(
                $"[{label}] scans={nScans} env={totalEnv} traces={traces.Count} features={features.Count} multi={multi} | " +
                $"BuildTraces={buildMs:F0}ms  Correct={corrMs:F0}ms  BuildFeatures={featMs:F0}ms");

            // Generous absolute tripwire: a representative TD size should be far under this.
            Assert.That(buildMs, Is.LessThan(60000), "BuildTraces exceeded 60 s tripwire");

            return new StageResult(totalEnv, traces.Count, features.Count, multi, buildMs, corrMs, featMs);
        }

        /// <summary>
        /// Synthetic generator: a rolling population of distinct neutral-mass species,
        /// each eluting over a contiguous scan window at several charges, with occasional
        /// +1 off-by-one twins. Seeded for reproducibility. <paramref name="activeSpeciesPerScan"/>
        /// controls how many species are live at once -- the open-trace count that drives
        /// BuildTraces cost -- and is held constant as nScans grows so scaling is meaningful.
        /// </summary>
        private static (MsDataScan[] scans, List<IReadOnlyList<IsotopicEnvelope>> perScan, int totalEnv)
            Generate(int nScans, int activeSpeciesPerScan, int chargesPerSpecies)
        {
            var rng = new Random(12345);
            const int elutionWidth = 20;

            // Size the species pool so ~activeSpeciesPerScan are live in any given scan.
            int nSpecies = (int)Math.Ceiling((double)nScans / elutionWidth) * activeSpeciesPerScan + activeSpeciesPerScan;
            var startScan = new int[nSpecies];
            var baseMass = new double[nSpecies];
            var baseCharge = new int[nSpecies];
            for (int s = 0; s < nSpecies; s++)
            {
                startScan[s] = rng.Next(0, Math.Max(1, nScans));
                baseMass[s] = 2000.0 + rng.NextDouble() * 28000.0; // 2k–30k Da
                baseCharge[s] = 6 + rng.Next(0, 12);               // lowest charge ~6..17
            }

            var scans = new MsDataScan[nScans];
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>(nScans);
            int totalEnv = 0;

            for (int i = 0; i < nScans; i++)
            {
                scans[i] = Scan(i + 1, i * 0.02); // ~0.02 min/scan
                var env = new List<IsotopicEnvelope>();
                for (int s = 0; s < nSpecies; s++)
                {
                    if (i < startScan[s] || i >= startScan[s] + elutionWidth) continue;
                    double intensity = 1e6 * (1 + rng.NextDouble());
                    for (int c = 0; c < chargesPerSpecies; c++)
                        env.Add(new IsotopicEnvelope(baseMass[s], intensity, baseCharge[s] + c));
                    if (rng.NextDouble() < 0.1) // occasional +1 twin at the lowest charge
                        env.Add(new IsotopicEnvelope(baseMass[s] + 1.00335, intensity * 0.5, baseCharge[s]));
                }
                totalEnv += env.Count;
                perScan.Add(env);
            }
            return (scans, perScan, totalEnv);
        }

        private static MsDataScan Scan(int oneBased, double rt)
            => new MsDataScan(
                new MzSpectrum(new[] { 500.0 }, new[] { 1.0 }, false),
                oneBased, msnOrder: 1, isCentroid: false, Polarity.Positive,
                retentionTime: rt, scanWindowRange: new MzRange(300, 2000),
                scanFilter: "ms1", mzAnalyzer: MZAnalyzerType.Unknown,
                totalIonCurrent: 1.0, injectionTime: null, noiseData: null, nativeId: null);
    }
}
