using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using MassSpectrometry.Deconvolution.FeatureTracing;
using MzLibUtil;
using NUnit.Framework;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Tests the pluggable feature-tracer boundary (<see cref="IMassFeatureTracer"/>) and the
    /// consensus-pipeline adapter (<see cref="ConsensusMassFeatureTracer"/>). The adapter must
    /// reproduce the direct MassTraceBuilder → TraceCorrector → MassFeatureBuilder recipe
    /// exactly, so MetaFlashDecon and consensus tracers are interchangeable on identical input.
    /// </summary>
    [TestFixture]
    public class TestMassFeatureTracer
    {
        private const double TolDa = 1.5;
        private const int MaxGap = 1;
        private const double CrossChargePpm = 10.0;

        [Test]
        public void ConsensusTracer_MatchesDirectPipeline()
        {
            var (scans, perScan) = BuildSample();

            var viaTracer = new ConsensusMassFeatureTracer(TolDa, MaxGap, CrossChargePpm)
                .TraceFeatures(scans, perScan);

            // The exact recipe the adapter wraps.
            var traces = MassTraceBuilder.BuildTraces(scans, perScan, TolDa, MaxGap);
            var corrected = traces.SelectMany(t => TraceCorrector.Correct(t)).ToList();
            var direct = MassFeatureBuilder.BuildFeatures(corrected, CrossChargePpm);

            Assert.That(viaTracer.Count, Is.EqualTo(direct.Count));
            Assert.That(viaTracer.Select(f => f.ConsensusMass).OrderBy(m => m).ToArray(),
                Is.EqualTo(direct.Select(f => f.ConsensusMass).OrderBy(m => m).ToArray()).Within(1e-9));
            Assert.That(viaTracer.Select(f => f.ChargeCount).OrderBy(c => c).ToArray(),
                Is.EqualTo(direct.Select(f => f.ChargeCount).OrderBy(c => c).ToArray()));
        }

        [Test]
        public void ConsensusTracer_StitchesChargesIntoOneFeature()
        {
            // Same neutral mass at two charges across two scans: charge-locked traces are
            // built per charge, then stitched cross-charge into ONE feature (ChargeCount 2).
            var (scans, perScan) = BuildSample();

            var features = new ConsensusMassFeatureTracer(TolDa, MaxGap, CrossChargePpm)
                .TraceFeatures(scans, perScan);

            Assert.That(features, Has.Count.EqualTo(1));
            Assert.That(features[0].ChargeCount, Is.EqualTo(2));
            Assert.That(features[0].ConsensusMass, Is.EqualTo(5000.0).Within(1e-6));
        }

        [Test]
        public void ConsensusTracer_IsAnIMassFeatureTracer()
        {
            IMassFeatureTracer tracer = new ConsensusMassFeatureTracer();
            Assert.That(tracer, Is.InstanceOf<IMassFeatureTracer>());
        }

        // ── MetaFlashDecon native tracer (neutral-mass tracing) ──────────────

        [Test]
        public void NativeTracer_CollapsesChargesIntoOneNeutralMassFeature()
        {
            // Same neutral mass at two charges, two scans: FLASHDeconv collapses charges per
            // scan, traces over neutral mass, then re-attaches charges -> one feature, ChargeCount 2.
            var (scans, perScan) = BuildSample();
            var features = new MetaFlashDeconMassFeatureTracer(
                    massTolerancePpm: 10, minTraceLengthRt: 0.1, minSampleRate: 0.05, maxOutlierScans: 1)
                .TraceFeatures(scans, perScan);

            Assert.That(features, Has.Count.EqualTo(1));
            Assert.That(features[0].ChargeCount, Is.EqualTo(2));
            Assert.That(features[0].ConsensusMass, Is.EqualTo(5000.0).Within(1e-6));
        }

        [Test]
        public void NativeTracer_TooShortTrace_IsFilteredOut()
        {
            // RT span is 0.2 (10.0 -> 10.2); a 1.0 minimum rejects it (FLASHDeconv length filter).
            var (scans, perScan) = BuildSample();
            var features = new MetaFlashDeconMassFeatureTracer(minTraceLengthRt: 1.0)
                .TraceFeatures(scans, perScan);
            Assert.That(features, Is.Empty);
        }

        [Test]
        public void NativeTracer_SeparatesDistinctMasses()
        {
            var (scans, perScan) = BuildTwoMassSample();
            var features = new MetaFlashDeconMassFeatureTracer(minTraceLengthRt: 0.1, maxOutlierScans: 1)
                .TraceFeatures(scans, perScan);

            Assert.That(features, Has.Count.EqualTo(2));
            Assert.That(features.Select(f => f.ConsensusMass).OrderBy(m => m).ToArray(),
                Is.EqualTo(new[] { 5000.0, 8000.0 }).Within(1e-6));
        }

        [Test]
        public void NativeTracer_IsAnIMassFeatureTracer()
        {
            IMassFeatureTracer tracer = new MetaFlashDeconMassFeatureTracer();
            Assert.That(tracer, Is.InstanceOf<IMassFeatureTracer>());
        }

        // Same neutral mass (5000 Da) at charges 5 and 6, present in two consecutive MS1 scans.
        private static (IReadOnlyList<MsDataScan>, IReadOnlyList<IReadOnlyList<IsotopicEnvelope>>) BuildSample()
        {
            var scans = new[] { Scan(1, 10.0), Scan(2, 10.2) };
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>
            {
                new[] { Env(5000.0, 1e7, 5), Env(5000.0, 1e7, 6) },
                new[] { Env(5000.0, 1e7, 5), Env(5000.0, 1e7, 6) },
            };
            return (scans, perScan);
        }

        // Two well-separated neutral masses (5000 @ z5, 8000 @ z6) in two consecutive scans.
        private static (IReadOnlyList<MsDataScan>, IReadOnlyList<IReadOnlyList<IsotopicEnvelope>>) BuildTwoMassSample()
        {
            var scans = new[] { Scan(1, 10.0), Scan(2, 10.2) };
            var perScan = new List<IReadOnlyList<IsotopicEnvelope>>
            {
                new[] { Env(5000.0, 1e7, 5), Env(8000.0, 1e7, 6) },
                new[] { Env(5000.0, 1e7, 5), Env(8000.0, 1e7, 6) },
            };
            return (scans, perScan);
        }

        private static IsotopicEnvelope Env(double mono, double intensity, int charge)
            => new IsotopicEnvelope(mono, intensity, charge);

        private static MsDataScan Scan(int oneBased, double rt)
            => new MsDataScan(
                new MzSpectrum(new[] { 500.0 }, new[] { 1.0 }, false),
                oneBased, msnOrder: 1, isCentroid: false, Polarity.Positive,
                retentionTime: rt, scanWindowRange: new MzRange(300, 2000),
                scanFilter: "ms1", mzAnalyzer: MZAnalyzerType.Unknown,
                totalIonCurrent: 1.0, injectionTime: null, noiseData: null, nativeId: null);
    }
}
