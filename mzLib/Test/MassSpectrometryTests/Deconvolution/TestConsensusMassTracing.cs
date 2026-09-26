using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MzLibUtil;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Step-1 SEMANTICS PROBES for the consensus mass-tracing engines
    /// (plan: E:\CodeReview\PR1069_test_and_perf_plan.md). These are confidence
    /// checks on the intended behaviour BEFORE the full lock-in suite -- if one
    /// fails we have found a bug to fix rather than a behaviour to cement.
    ///
    /// Probes 1-3 assert behaviour we are confident is correct (charge-lock,
    /// off-by-one rescue, cross-charge stitching). Probe 4 is OBSERVATIONAL: it
    /// characterises a known-open design question (the loose 1.5 Da Pass-B window
    /// merging near-isobaric same-charge species) for domain review -- the
    /// behaviour here is the question, not the locked answer.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestConsensusMassTracing
    {
        // Pass-B (loose) grouping params used by the consensus pipeline.
        private const double PassBToleranceDa = 1.5;
        private const int MaxGap = 1;
        private const double CrossChargeMassPpm = 10.0;

        [Test]
        public void Probe_BuildTraces_IsChargeLocked()
        {
            // Same neutral mass at two charges across two scans. Charge-locked
            // grouping must yield ONE trace per charge, never a cross-charge merge.
            var scans = new[] { Scan(1, 10.0), Scan(2, 10.2) };
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>
            {
                new[] { Env(5000.0, 1e7, 5), Env(5000.0, 1e7, 6) },
                new[] { Env(5000.0, 1e7, 5), Env(5000.0, 1e7, 6) },
            };

            var traces = MassTraceBuilder.BuildTraces(scans, perScan, PassBToleranceDa, MaxGap);

            Assert.That(traces.Count, Is.EqualTo(2), "one trace per charge");
            Assert.That(traces.Select(t => t.Charge).OrderBy(c => c), Is.EqualTo(new[] { 5, 6 }));
            Assert.That(traces.All(t => t.Envelopes.Count == 2), "each charge seen in both scans");
        }

        [Test]
        public void Probe_OffByOne_IsRescuedToConsensus()
        {
            // Three envelopes at the true mono mass plus one +1-isotope shadow, same
            // charge, in separate scans so the loose grouper cohorts all four into one
            // trace (the one-per-scan guard forbids two envelopes joining in one scan).
            const double mono = 5000.0;
            double shadow = mono + TraceCorrector.IsotopeSpacingDa; // +1.00335 Da
            var scans = new[] { Scan(1, 10.0), Scan(2, 10.2), Scan(3, 10.4), Scan(4, 10.6) };
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>
            {
                new[] { Env(mono, 1e7, 10) },
                new[] { Env(mono, 1e7, 10) },
                new[] { Env(mono, 1e7, 10) },
                new[] { Env(shadow, 1e7, 10) },
            };

            var traces = MassTraceBuilder.BuildTraces(scans, perScan, PassBToleranceDa, MaxGap);
            Assert.That(traces.Count, Is.EqualTo(1), "loose grouping cohorts the +1 twin with its parent");

            // A +1.00335 shadow is an integer isotope offset, so it is merged+corrected,
            // never split -- Single() also asserts that no spurious split happened.
            var corrected = TraceCorrector.Correct(traces[0]).Single();

            Assert.That(corrected.ConsensusMass, Is.EqualTo(mono).Within(1e-6), "median anchors on the true mono");
            Assert.That(corrected.CorrectionCount, Is.EqualTo(1), "exactly the shadow envelope corrected");
            Assert.That(corrected.CorrectedSpread, Is.LessThan(corrected.OriginalSpread), "spread shrinks after rescue");

            var fixedEnv = corrected.Envelopes.Single(e => e.WasCorrected);
            Assert.That(fixedEnv.OriginalMass, Is.EqualTo(shadow).Within(1e-6));
            Assert.That(fixedEnv.CorrectedMass, Is.EqualTo(mono).Within(1e-6));
        }

        [Test]
        public void Probe_CrossCharge_StitchesSameMassDifferentCharge()
        {
            // One species seen at z5 and z6 (identical neutral mass), overlapping RT.
            var t5 = MakeCorrectedTrace(charge: 5, consensus: 8000.0, rtStart: 20.0, rtEnd: 21.0, intensity: 2e7);
            var t6 = MakeCorrectedTrace(charge: 6, consensus: 8000.0, rtStart: 20.1, rtEnd: 21.1, intensity: 1.5e7);

            var features = MassFeatureBuilder.BuildFeatures(new[] { t5, t6 }, CrossChargeMassPpm);

            Assert.That(features.Count, Is.EqualTo(1), "same mass + RT overlap + differing charge => one feature");
            Assert.That(features[0].ChargeCount, Is.EqualTo(2));
            Assert.That(features[0].Charges.OrderBy(c => c), Is.EqualTo(new[] { 5, 6 }));
        }

        [Test]
        public void Probe_NearIsobaricSameCharge_OPENQUESTION_CharacteriseFalseMerge()
        {
            // OPEN DESIGN QUESTION (domain review): two DISTINCT species 0.8 Da apart at
            // the SAME charge. When their scan presence interleaves (B appears in a scan
            // where A is absent), the loose 1.5 Da anchor window lets B join A's trace --
            // a false merge that swallows the 0.8-Da-distinct species. This probe
            // documents the GROUPING stage: the loose grouper intentionally cohorts the
            // two species so TraceCorrector can decide. DECIDED (option 2): TraceCorrector
            // splits resolvable distinct species back out -- see Probe_Deamidation_HighRes_*.
            var scans = new[] { Scan(1, 10.0), Scan(2, 10.2), Scan(3, 10.4) };
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>
            {
                new[] { Env(5000.0, 1e7, 10) },   // species A
                new[] { Env(5000.8, 1e7, 10) },   // species B, A absent this scan
                new[] { Env(5000.0, 1e7, 10) },   // species A again
            };

            var traces = MassTraceBuilder.BuildTraces(scans, perScan, PassBToleranceDa, MaxGap);

            TestContext.Out.WriteLine($"[near-isobaric probe] traces={traces.Count}");
            foreach (var t in traces)
                TestContext.Out.WriteLine(
                    $"  charge={t.Charge} n={t.Envelopes.Count} spread={t.MassSpread:F4} " +
                    $"masses=[{string.Join(",", t.Envelopes.Select(e => e.Mass.ToString("F3")))}]");

            // Documents CURRENT behaviour: the interleaved B merges into A's trace.
            Assert.That(traces.Count, Is.EqualTo(1), "DOCUMENTS current false-merge; revisit per domain decision");
            Assert.That(traces[0].MassSpread, Is.EqualTo(0.8).Within(1e-6));
        }

        [Test]
        public void Probe_Deamidation_HighRes_SurfacesAsSeparateTrace()
        {
            // High resolution (exact masses): an unmodified cluster plus a deamidated
            // cluster at +0.98402 Da, same charge. 0.98402 differs from the 1.00335 Da
            // isotope spacing by 0.0193 Da; with tight scatter the adaptive window resolves
            // that gap, so the deamidated form SPLITS into its own corrected trace instead
            // of being snapped to the unmodified mass.
            const double mono = 8000.0;
            const double deam = mono + 0.98402;
            var trace = MakeMassTrace(charge: 8,
                masses: new[] { mono, mono, mono, mono, mono, deam, deam, deam });

            var corrected = TraceCorrector.Correct(trace);

            Assert.That(corrected.Count, Is.EqualTo(2), "deamidated species surfaces as its own trace");
            var unmod = corrected.Single(c => c.ConsensusMass < mono + 0.5);
            var modified = corrected.Single(c => c.ConsensusMass > mono + 0.5);
            Assert.That(unmod.ConsensusMass, Is.EqualTo(mono).Within(1e-6));
            Assert.That(modified.ConsensusMass, Is.EqualTo(deam).Within(1e-6));
            Assert.That(corrected.SelectMany(c => c.Envelopes).Any(e => e.WasCorrected), Is.False,
                "neither cluster is an isotope error, so nothing is mass-snapped");
        }

        [Test]
        public void Probe_Deamidation_LowRes_FallsBackToMerge()
        {
            // Same unmodified + deamidated split, but with ~0.03 Da mass scatter so the
            // estimated 3-sigma window exceeds the 0.0193 Da gap. The two can no longer be
            // told apart, so the conservative fallback treats the +~0.98 cluster as an
            // isotope error and merges it -- documenting that low-resolution data cannot
            // preserve deamidation (a resolution limit, not a defect).
            const double mono = 8000.0;
            var trace = MakeMassTrace(charge: 8, masses: new[]
            {
                mono - 0.03, mono + 0.02, mono - 0.01, mono + 0.03, mono,   // scattered base
                mono + 0.95, mono + 1.00, mono + 0.97,                      // scattered ~+0.98
            });

            var corrected = TraceCorrector.Correct(trace);

            Assert.That(corrected.Count, Is.EqualTo(1), "low-res cannot resolve deamidation from an isotope error");
        }

        // ──────────── step-3 fixes: determinism + input guards ────────────

        [Test]
        public void BuildFeatures_IsDeterministic_AcrossInputOrder()
        {
            // Same traces in two different input orders must produce identical features
            // (Id, consensus mass, charges). Guards the nondeterministic Dictionary
            // enumeration order that previously leaked into feature IDs and file row order.
            var traces = new[]
            {
                MakeCorrectedTrace(charge: 5, consensus: 8000.0, rtStart: 20.0, rtEnd: 21.0, intensity: 2e7),
                MakeCorrectedTrace(charge: 6, consensus: 8000.0, rtStart: 20.0, rtEnd: 21.0, intensity: 2e7),
                MakeCorrectedTrace(charge: 5, consensus: 3000.0, rtStart: 10.0, rtEnd: 11.0, intensity: 1e7),
                MakeCorrectedTrace(charge: 4, consensus: 12000.0, rtStart: 30.0, rtEnd: 31.0, intensity: 3e7),
            };

            var a = MassFeatureBuilder.BuildFeatures(traces, CrossChargeMassPpm);
            var b = MassFeatureBuilder.BuildFeatures(traces.Reverse().ToArray(), CrossChargeMassPpm);

            Assert.That(a.Count, Is.EqualTo(b.Count));
            for (int i = 0; i < a.Count; i++)
            {
                Assert.That(a[i].Id, Is.EqualTo(b[i].Id), $"feature {i} Id");
                Assert.That(a[i].ConsensusMass, Is.EqualTo(b[i].ConsensusMass).Within(1e-9), $"feature {i} mass");
                Assert.That(a[i].Charges.OrderBy(c => c), Is.EqualTo(b[i].Charges.OrderBy(c => c)), $"feature {i} charges");
            }
            Assert.That(a.Select(f => f.ConsensusMass), Is.Ordered, "features ordered by ascending consensus mass");
        }

        [Test]
        public void Correct_EmptyTrace_Throws()
            => Assert.Throws<System.ArgumentException>(
                () => TraceCorrector.Correct(new MassTrace { Id = 1, Charge = 5 }));

        [Test]
        public void Finalise_NoTraces_Throws()
            => Assert.Throws<System.InvalidOperationException>(() => new MassFeature().Finalise());

        [Test]
        public void BuildTraces_LengthMismatch_Throws()
        {
            var scans = new[] { Scan(1, 10.0), Scan(2, 10.2) };
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>> { new[] { Env(5000.0, 1e7, 5) } };
            Assert.Throws<System.ArgumentException>(
                () => MassTraceBuilder.BuildTraces(scans, perScan, PassBToleranceDa, MaxGap));
        }

        // ──────────── step-4 coverage: MassTraceBuilder ────────────

        [Test]
        public void BuildTraces_GapBoundary_RetiresWhenGapExceedsMaxGap()
        {
            // maxGap=1: a one-scan gap keeps the trace open; a two-scan gap retires it
            // and the next envelope starts a fresh trace.
            const double m = 5000.0; const int z = 10;
            var withinGap = MassTraceBuilder.BuildTraces(
                new[] { Scan(1, 1.0), Scan(2, 1.1), Scan(3, 1.2) },
                new List<IReadOnlyList<IsotopicEnvelope>>
                {
                    new[] { Env(m, 1e7, z) },
                    System.Array.Empty<IsotopicEnvelope>(),                 // 1-scan gap
                    new[] { Env(m, 1e7, z) },
                },
                PassBToleranceDa, MaxGap);
            Assert.That(withinGap.Count, Is.EqualTo(1), "gap == maxGap keeps the trace open");

            var beyondGap = MassTraceBuilder.BuildTraces(
                new[] { Scan(1, 1.0), Scan(2, 1.1), Scan(3, 1.2), Scan(4, 1.3) },
                new List<IReadOnlyList<IsotopicEnvelope>>
                {
                    new[] { Env(m, 1e7, z) },
                    System.Array.Empty<IsotopicEnvelope>(),
                    System.Array.Empty<IsotopicEnvelope>(),                 // 2-scan gap
                    new[] { Env(m, 1e7, z) },
                },
                PassBToleranceDa, MaxGap);
            Assert.That(beyondGap.Count, Is.EqualTo(2), "gap > maxGap retires and starts a new trace");
        }

        [Test]
        public void BuildTraces_BestMatchWins_JoinsClosestAnchor()
        {
            // Two open same-charge traces (anchors 2 Da apart). An envelope within tolerance
            // of both joins the closer anchor.
            var traces = MassTraceBuilder.BuildTraces(
                new[] { Scan(1, 1.0), Scan(2, 1.1) },
                new List<IReadOnlyList<IsotopicEnvelope>>
                {
                    new[] { Env(5000.0, 1e7, 10), Env(5002.0, 1e7, 10) }, // two separate traces
                    new[] { Env(5000.8, 1e7, 10) },                        // closer to 5000.0
                },
                PassBToleranceDa, MaxGap);

            Assert.That(traces.Count, Is.EqualTo(2));
            var joined = traces.Single(t => t.Envelopes.Count == 2);
            Assert.That(joined.AnchorMass, Is.EqualTo(5000.0).Within(1e-9), "joined the closer anchor");
        }

        [Test]
        public void BuildTraces_TwoEnvelopesSameScan_DoNotJoinSameTrace()
        {
            var traces = MassTraceBuilder.BuildTraces(
                new[] { Scan(1, 1.0) },
                new List<IReadOnlyList<IsotopicEnvelope>>
                {
                    new[] { Env(5000.0, 1e7, 10), Env(5000.1, 1e7, 10) },
                },
                PassBToleranceDa, MaxGap);
            Assert.That(traces.Count, Is.EqualTo(2), "two same-scan envelopes cannot join one trace");
            Assert.That(traces.All(t => t.Envelopes.Count == 1));
        }

        [Test]
        public void BuildTraces_AnchorFixed_DriftBeyondToleranceStartsNewTrace()
        {
            // Anchor stays at the first mass. Masses drifting +0.6/scan stay in the trace
            // until one lands beyond anchor +/- tolerance, which starts a new trace.
            var traces = MassTraceBuilder.BuildTraces(
                new[] { Scan(1, 1.0), Scan(2, 1.1), Scan(3, 1.2), Scan(4, 1.3) },
                new List<IReadOnlyList<IsotopicEnvelope>>
                {
                    new[] { Env(5000.0, 1e7, 10) },
                    new[] { Env(5000.6, 1e7, 10) },
                    new[] { Env(5001.2, 1e7, 10) },
                    new[] { Env(5001.8, 1e7, 10) }, // 1.8 from the fixed anchor (5000) -> falls out
                },
                PassBToleranceDa, MaxGap);
            Assert.That(traces.Count, Is.EqualTo(2), "drift beyond tolerance of the fixed anchor starts a new trace");
            Assert.That(traces.Max(t => t.Envelopes.Count), Is.EqualTo(3));
        }

        [Test]
        public void BuildTraces_NoScans_ReturnsEmpty()
            => Assert.That(
                MassTraceBuilder.BuildTraces(
                    System.Array.Empty<MsDataScan>(),
                    new List<IReadOnlyList<IsotopicEnvelope>>(),
                    PassBToleranceDa, MaxGap),
                Is.Empty);

        // ──────────── step-4 coverage: TraceCorrector ────────────

        [Test]
        public void Correct_MinusOneShadow_IsRescued()
        {
            const double mono = 5000.0;
            var trace = MakeMassTrace(charge: 10,
                masses: new[] { mono, mono, mono, mono - TraceCorrector.IsotopeSpacingDa });
            var corrected = TraceCorrector.Correct(trace).Single();

            Assert.That(corrected.ConsensusMass, Is.EqualTo(mono).Within(1e-6));
            Assert.That(corrected.CorrectionCount, Is.EqualTo(1));
            var fixedEnv = corrected.Envelopes.Single(e => e.WasCorrected);
            Assert.That(fixedEnv.OriginalMass, Is.EqualTo(mono - TraceCorrector.IsotopeSpacingDa).Within(1e-6));
            Assert.That(fixedEnv.CorrectedMass, Is.EqualTo(mono).Within(1e-6));
        }

        [Test]
        public void Correct_AllZeroWeights_FallsBackToUnweightedMedian()
        {
            // IntensityWeight over zero-intensity envelopes => total weight 0 => the weighted
            // median degenerates to the plain median (exercises the degenerate branch and the
            // IntensityWeight delegate).
            var trace = MakeMassTrace(charge: 10, masses: new[] { 5000.0, 5000.0, 5000.0 }, intensity: 0.0);
            var corrected = TraceCorrector.Correct(trace, TraceCorrector.IntensityWeight);
            Assert.That(corrected.Count, Is.EqualTo(1));
            Assert.That(corrected[0].ConsensusMass, Is.EqualTo(5000.0).Within(1e-9));
        }

        // ──────────── step-4 coverage: MassFeatureBuilder ────────────

        [Test]
        public void BuildFeatures_PpmWindow_MergesInsideSeparatesOutside()
        {
            // 10 ppm at 8000 Da ~ 0.08 Da window.
            var inside = MassFeatureBuilder.BuildFeatures(new[]
            {
                MakeCorrectedTrace(5, 8000.00, 20.0, 21.0, 1e7),
                MakeCorrectedTrace(6, 8000.05, 20.0, 21.0, 1e7),
            }, CrossChargeMassPpm);
            Assert.That(inside.Count, Is.EqualTo(1), "0.05 Da apart (< 0.08) merges");

            var outside = MassFeatureBuilder.BuildFeatures(new[]
            {
                MakeCorrectedTrace(5, 8000.00, 20.0, 21.0, 1e7),
                MakeCorrectedTrace(6, 8000.10, 20.0, 21.0, 1e7),
            }, CrossChargeMassPpm);
            Assert.That(outside.Count, Is.EqualTo(2), "0.10 Da apart (> 0.08) stays separate");
        }

        [Test]
        public void BuildFeatures_SameCharge_NotMerged()
        {
            var features = MassFeatureBuilder.BuildFeatures(new[]
            {
                MakeCorrectedTrace(5, 8000.0, 20.0, 21.0, 1e7),
                MakeCorrectedTrace(5, 8000.0, 20.0, 21.0, 1e7),
            }, CrossChargeMassPpm);
            Assert.That(features.Count, Is.EqualTo(2), "same-charge co-mass traces are not stitched");
        }

        [Test]
        public void BuildFeatures_RtDisjoint_NotMerged()
        {
            var features = MassFeatureBuilder.BuildFeatures(new[]
            {
                MakeCorrectedTrace(5, 8000.0, 10.0, 11.0, 1e7),
                MakeCorrectedTrace(6, 8000.0, 20.0, 21.0, 1e7),
            }, CrossChargeMassPpm);
            Assert.That(features.Count, Is.EqualTo(2), "mass agrees but RT windows are disjoint -> not stitched");
        }

        [Test]
        public void BuildFeatures_TransitiveUnion_ChainsIntoOneFeature()
        {
            // A~B and B~C are each within the ppm window, but A and C are beyond it:
            // union-find must still place all three in one feature.
            var features = MassFeatureBuilder.BuildFeatures(new[]
            {
                MakeCorrectedTrace(5, 8000.00, 20.0, 21.0, 1e7),
                MakeCorrectedTrace(6, 8000.07, 20.0, 21.0, 1e7),
                MakeCorrectedTrace(7, 8000.14, 20.0, 21.0, 1e7),
            }, CrossChargeMassPpm);
            Assert.That(features.Count, Is.EqualTo(1));
            Assert.That(features[0].Charges.OrderBy(c => c), Is.EqualTo(new[] { 5, 6, 7 }));
        }

        [Test]
        public void BuildFeatures_NoTraces_ReturnsEmpty()
            => Assert.That(
                MassFeatureBuilder.BuildFeatures(System.Array.Empty<CorrectedTrace>(), CrossChargeMassPpm),
                Is.Empty);

        // ──────────── step-4 coverage: MassFeature ────────────

        [Test]
        public void Finalise_ZeroIntensityTraces_UsesPlainAverageConsensus()
        {
            // Every constituent trace has zero total intensity, so the intensity-weighted
            // mean is undefined and Finalise falls back to a plain average of consensus masses.
            var feature = new MassFeature
            {
                Traces =
                {
                    MakeCorrectedTrace(5, 1000.0, 10.0, 11.0, 0.0),
                    MakeCorrectedTrace(6, 2000.0, 10.0, 11.0, 0.0),
                },
            };
            feature.Finalise();
            Assert.That(feature.ConsensusMass, Is.EqualTo(1500.0).Within(1e-9));
        }

        // ───────────────────────── helpers ─────────────────────────

        private static IsotopicEnvelope Env(double mono, double intensity, int charge)
            => new IsotopicEnvelope(mono, intensity, charge);

        private static MsDataScan Scan(int oneBased, double rt)
            => new MsDataScan(
                new MzSpectrum(new[] { 500.0 }, new[] { 1.0 }, false),
                oneBased, msnOrder: 1, isCentroid: false, Polarity.Positive,
                retentionTime: rt, scanWindowRange: new MzRange(300, 2000),
                scanFilter: "ms1", mzAnalyzer: MZAnalyzerType.Unknown,
                totalIonCurrent: 1.0, injectionTime: null, noiseData: null, nativeId: null);

        private static MassTrace MakeMassTrace(int charge, double[] masses, double intensity = 1e7)
        {
            var mt = new MassTrace { Id = 1, Charge = charge, AnchorMass = masses[0] };
            for (int i = 0; i < masses.Length; i++)
                mt.Envelopes.Add((i, i + 1, 10.0 + i * 0.1, masses[i], intensity));
            return mt;
        }

        private static CorrectedTrace MakeCorrectedTrace(int charge, double consensus, double rtStart, double rtEnd, double intensity)
        {
            var ct = new CorrectedTrace { Id = charge, Charge = charge, ConsensusMass = consensus };
            ct.Envelopes.Add(new CorrectedEnvelope
            {
                ScanIndex = 0, ScanNumber = 1, RT = rtStart, OriginalMass = consensus,
                CorrectedMass = consensus, Charge = charge, Intensity = intensity / 2,
            });
            ct.Envelopes.Add(new CorrectedEnvelope
            {
                ScanIndex = 1, ScanNumber = 2, RT = rtEnd, OriginalMass = consensus,
                CorrectedMass = consensus, Charge = charge, Intensity = intensity / 2,
            });
            return ct;
        }
    }
}
