using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// HEAD-TO-HEAD: the PR's new consensus tracer (<see cref="MassTraceBuilder.BuildTraces"/>)
    /// vs. mzLib's existing XIC machinery (<see cref="MassIndexingEngine"/> +
    /// <see cref="IndexingEngine{T}.GetAllXics"/>) -- the reuse question nic raised on PR #1069.
    ///
    /// Both consume the SAME deconvoluted per-scan envelopes (decon is run once, up front), so
    /// the timings isolate the TRACING algorithm, not deconvolution. The existing engine is fed
    /// the identical envelopes through a test-only injection seam (no production change), and is
    /// given an <see cref="AbsoluteTolerance"/> of 1.5 Da + missed-scans = MaxGap so its tracing
    /// window matches the consensus tracer's (which is deliberately wide so TraceCorrector can
    /// later sort out off-by-one isotope errors).
    ///
    /// Reports, for each engine: wall-clock, bytes allocated, trace count, trace-length
    /// distribution, and a grouping-agreement metric (how aligned the two groupings are).
    /// This answers PERFORMANCE + grouping similarity. End-to-end PROTEOFORM YIELD (the
    /// decisive q&lt;0.01 metric) is a separate pipeline run, not this micro-benchmark.
    ///
    /// Run: dotnet test --filter FullyQualifiedName~TestXicVsMassTraceComparison
    /// </summary>
    [TestFixture]
    [Explicit("Head-to-head perf/grouping comparison on a local Jurkat TD mzML; needs E:\\Projects\\JurkatTopDown")]
    [ExcludeFromCodeCoverage]
    public class TestXicVsMassTraceComparison
    {
        // One representative top-down file. Override via the env var to sweep others.
        private static string MzmlPath =>
            Environment.GetEnvironmentVariable("XIC_CMP_MZML")
            ?? @"E:\Projects\JurkatTopDown\02-17-20_jurkat_td_rep1_fract1.mzML";

        // Classic top-down decon parameters (match JurkatMs1FeatureGenerator).
        private const int MinCharge = 1;
        private const int MaxCharge = 60;
        private const double DeconPpm = 4.0;
        private const double IntensityRatio = 3.0;

        // Consensus-tracing parameters (match the end-to-end driver).
        private const double TraceToleranceDa = 1.5;
        private const int MaxGap = 1;

        [Test]
        public void CompareTracers_OnRealTopDownFile()
        {
            Assume.That(System.IO.File.Exists(MzmlPath), $"mzML not found: {MzmlPath}");

            // ---- shared input: deconvolute every MS1 scan once ----
            var dataFile = MsDataFileReader.GetDataFile(MzmlPath).LoadAllStaticData();
            var ms1Scans = dataFile.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToList();

            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: MinCharge, maxCharge: MaxCharge,
                deconPpm: DeconPpm, intensityRatio: IntensityRatio);

            var perScanEnvelopes = ms1Scans
                .Select(s => (IReadOnlyList<IsotopicEnvelope>)Deconvoluter
                    .Deconvolute(s.MassSpectrum, deconParams).ToList())
                .ToList();
            int totalEnv = perScanEnvelopes.Sum(s => s.Count);
            double maxMass = perScanEnvelopes.SelectMany(s => s).DefaultIfEmpty().Max(e => e?.MonoisotopicMass ?? 0);

            TestContext.Progress.WriteLine(
                $"[input] {System.IO.Path.GetFileName(MzmlPath)}  ms1Scans={ms1Scans.Count}  envelopes={totalEnv}  maxMass={maxMass:F0} Da");

            // ---- warm up both paths (JIT) on the real data, results discarded ----
            _ = MassTraceBuilder.BuildTraces(ms1Scans, perScanEnvelopes, TraceToleranceDa, MaxGap);
            var warmEngine = BuildInjectedEngine(ms1Scans, perScanEnvelopes, maxMass);
            _ = warmEngine.GetAllXics(new AbsoluteTolerance(TraceToleranceDa), MaxGap, double.MaxValue, numPeakThreshold: 1);

            // ---- Path A: consensus MassTraceBuilder ----
            GcReset();
            long aBytes0 = GC.GetTotalAllocatedBytes(true);
            var swA = Stopwatch.StartNew();
            var traces = MassTraceBuilder.BuildTraces(ms1Scans, perScanEnvelopes, TraceToleranceDa, MaxGap);
            swA.Stop();
            long aBytes = GC.GetTotalAllocatedBytes(true) - aBytes0;

            // ---- Path B: existing MassIndexingEngine + GetAllXics ----
            // Index-build (binning the same envelopes) is part of B's required cost; time it separately.
            GcReset();
            long bIdxBytes0 = GC.GetTotalAllocatedBytes(true);
            var swBidx = Stopwatch.StartNew();
            var engine = BuildInjectedEngine(ms1Scans, perScanEnvelopes, maxMass);
            swBidx.Stop();
            long bIdxBytes = GC.GetTotalAllocatedBytes(true) - bIdxBytes0;

            GcReset();
            long bBytes0 = GC.GetTotalAllocatedBytes(true);
            var swB = Stopwatch.StartNew();
            var xics = engine.GetAllXics(new AbsoluteTolerance(TraceToleranceDa), MaxGap, double.MaxValue, numPeakThreshold: 1);
            swB.Stop();
            long bBytes = GC.GetTotalAllocatedBytes(true) - bBytes0;

            // ---- shape metrics ----
            var aLens = traces.Select(t => t.Envelopes.Count).ToList();
            var bLens = xics.Select(x => x.Peaks.Count).ToList();
            int aClaimed = aLens.Sum();
            int bClaimed = bLens.Sum();

            // ---- grouping agreement ----
            // Key each envelope occurrence by (scanIndex, charge, mass to 3dp). Label it with the
            // trace id assigned by each engine, then count distinct (idA,idB) pairs: a perfectly
            // aligned grouping yields ~one pair per trace, so distinctPairs near max(traceCounts)
            // means high agreement; far above means the engines split/merge differently.
            var aLabel = new Dictionary<(int, int, long), int>(totalEnv);
            for (int t = 0; t < traces.Count; t++)
                foreach (var e in traces[t].Envelopes)
                    aLabel[Key(e.ScanIndex, traces[t].Charge, e.Mass)] = t;

            var bLabel = new Dictionary<(int, int, long), int>(totalEnv);
            for (int x = 0; x < xics.Count; x++)
                foreach (var p in xics[x].Peaks)
                    bLabel[Key(p.ZeroBasedScanIndex, (p as IndexedMass)?.Charge ?? 0, p.M)] = x;

            int shared = 0, distinctPairs;
            var pairs = new HashSet<(int, int)>();
            foreach (var kv in aLabel)
                if (bLabel.TryGetValue(kv.Key, out int bId)) { shared++; pairs.Add((kv.Value, bId)); }
            distinctPairs = pairs.Count;

            // ---- report ----
            void Report(string label, double ms, long bytes, int nTraces, List<int> lens, int claimed)
            {
                double mean = lens.Count == 0 ? 0 : lens.Average();
                int singles = lens.Count(l => l == 1);
                int multi = lens.Count(l => l >= 2);
                int maxLen = lens.Count == 0 ? 0 : lens.Max();
                TestContext.Out.WriteLine(
                    $"  {label,-26} time={ms,8:F1} ms  alloc={bytes / 1_048_576.0,7:F1} MB  " +
                    $"traces={nTraces,7}  claimed={claimed,7}  len(mean={mean:F2} max={maxLen} singles={singles} multi={multi})");
            }

            TestContext.Out.WriteLine($"\n=== XIC vs MassTrace on {System.IO.Path.GetFileName(MzmlPath)} ===");
            TestContext.Out.WriteLine($"  ms1Scans={ms1Scans.Count}  envelopes={totalEnv}  maxMass={maxMass:F0} Da  (tol={TraceToleranceDa} Da, maxGap={MaxGap})");
            Report("A: MassTraceBuilder", swA.Elapsed.TotalMilliseconds, aBytes, traces.Count, aLens, aClaimed);
            Report("B: GetAllXics", swB.Elapsed.TotalMilliseconds, bBytes, xics.Count, bLens, bClaimed);
            TestContext.Out.WriteLine(
                $"  B index-build (extra):     time={swBidx.Elapsed.TotalMilliseconds,8:F1} ms  alloc={bIdxBytes / 1_048_576.0,7:F1} MB");
            TestContext.Out.WriteLine(
                $"  B total (index+trace):     time={(swBidx.Elapsed.TotalMilliseconds + swB.Elapsed.TotalMilliseconds),8:F1} ms");
            TestContext.Out.WriteLine(
                $"  grouping agreement: sharedEnvelopes={shared}/{totalEnv}  distinct(idA,idB) pairs={distinctPairs}  " +
                $"(A traces={traces.Count}, B traces={xics.Count}; closer to max => more aligned)");

            // Sanity, not a perf gate: both engines must claim essentially every envelope (both keep singletons).
            Assert.That(aClaimed, Is.EqualTo(totalEnv), "MassTraceBuilder dropped envelopes");
            Assert.That(bClaimed, Is.GreaterThan(0), "GetAllXics produced no peaks");
        }

        private static (int, int, long) Key(int scanIndex, int charge, double mass)
            => (scanIndex, charge, (long)Math.Round(mass * 1000.0));

        private static void GcReset()
        {
            GC.Collect();
            GC.WaitForPendingFinalizers();
            GC.Collect();
        }

        /// <summary>
        /// Build a MassIndexingEngine populated from pre-deconvoluted envelopes (the same ones the
        /// consensus tracer sees), bypassing the engine's internal per-scan deconvolution so the
        /// comparison is decon-free. MaxMass is raised above the observed maximum so the engine's
        /// default 30 kDa bin cap is not the thing that differs.
        /// </summary>
        private static InjectedMassIndexingEngine BuildInjectedEngine(
            IReadOnlyList<MsDataScan> ms1Scans,
            IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes,
            double maxMass)
        {
            var engine = new InjectedMassIndexingEngine { MaxMass = (int)Math.Ceiling(maxMass) + 2 };
            engine.IndexFromEnvelopes(ms1Scans, perScanEnvelopes);
            return engine;
        }

        /// <summary>
        /// Test-only seam: exposes the protected index so we can inject shared envelopes. Mirrors
        /// <see cref="MassIndexingEngine.IndexPeaks"/> binning exactly (BinsPerDalton = 1), minus
        /// the deconvolution step.
        /// </summary>
        private sealed class InjectedMassIndexingEngine : MassIndexingEngine
        {
            public void IndexFromEnvelopes(
                IReadOnlyList<MsDataScan> ms1Scans,
                IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes)
            {
                IndexedPeaks = new List<IndexedMass>[MaxMass];
                ScanInfoArray = new ScanInfo[ms1Scans.Count];

                for (int scanIndex = 0; scanIndex < ms1Scans.Count; scanIndex++)
                {
                    var scan = ms1Scans[scanIndex];
                    ScanInfoArray[scanIndex] = new ScanInfo(scan.OneBasedScanNumber, scanIndex, scan.RetentionTime, scan.MsnOrder);

                    foreach (var envelope in perScanEnvelopes[scanIndex])
                    {
                        int roundedMass = (int)Math.Round(envelope.MonoisotopicMass, 0); // BinsPerDalton = 1
                        if (roundedMass < 0 || roundedMass >= IndexedPeaks.Length)
                            continue;
                        IndexedPeaks[roundedMass] ??= new List<IndexedMass>();
                        // Appended in ascending scan order, as GetXic's binary search requires.
                        IndexedPeaks[roundedMass].Add(new IndexedMass(envelope, scan.RetentionTime, scanIndex, scan.MsnOrder));
                    }
                }
            }
        }
    }
}
