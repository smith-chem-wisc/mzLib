using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Test.FileReadingTests;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.MassSpectrometryTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestFromFileDeconvolution
    {
        // -------- helpers --------------------------------------------------------

        private static MsDataScan MakeMs1(int oneBased, double rt) =>
            new MsDataScan(
                MinimalMs1Spectrum(),
                oneBased, msnOrder: 1, isCentroid: false, Polarity.Positive,
                retentionTime: rt, scanWindowRange: new MzRange(300, 2000),
                scanFilter: "ms1", mzAnalyzer: MZAnalyzerType.Unknown,
                totalIonCurrent: 1.0, injectionTime: null, noiseData: null, nativeId: null);

        private static MsDataScan MakeMs2(int oneBased, double rt, double? isolationMz, double? isolationWidth, int precursorScan = 1, int msnOrder = 2) =>
            new MsDataScan(
                MinimalMs2Spectrum(),
                oneBased, msnOrder: msnOrder, isCentroid: false, Polarity.Positive,
                retentionTime: rt, scanWindowRange: new MzRange(100, 2000),
                scanFilter: "ms2", mzAnalyzer: MZAnalyzerType.Unknown,
                totalIonCurrent: 1.0, injectionTime: null, noiseData: null, nativeId: null,
                selectedIonMz: isolationMz, selectedIonChargeStateGuess: null, selectedIonIntensity: null,
                isolationMZ: isolationMz, isolationWidth: isolationWidth,
                dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: precursorScan,
                selectedIonMonoisotopicGuessMz: null);

        private static MzSpectrum MinimalMs1Spectrum() => new MzSpectrum(new double[] { 500.0 }, new double[] { 1.0 }, false);
        private static MzSpectrum MinimalMs2Spectrum() => new MzSpectrum(new double[] { 100.0 }, new double[] { 1.0 }, false);

        // Feature centered at 600 m/z, charge 2 (mass ~ 1197.99), eluting from RT 10 to RT 15.
        private static SingleChargeMs1Feature Feature(double mz = 600.0, int charge = 2,
            double rtStart = 10.0, double rtEnd = 15.0, double intensity = 1e5)
            => new SingleChargeMs1Feature(mz, charge, rtStart, rtEnd, intensity);

        private static FromFileDeconvolutionParameters Params(params ISingleChargeMs1Feature[] feats)
            => new FromFileDeconvolutionParameters(feats, minCharge: 1, maxCharge: 60);

        // -------- direct algorithm tests via Deconvoluter.Deconvolute -----------
        // These tests bypass MsDataScan and exercise FromFileDeconvolutionAlgorithm
        // through Deconvoluter directly, with hand-built MzRtRange inputs.

        [Test]
        public void DirectAlgorithm_FeatureInsideRange_EmitsOneEnvelope()
        {
            var feat = Feature();
            var spectrum = MinimalMs1Spectrum();
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 12.0, maxRt: 13.0);

            var envelopes = Deconvoluter.Deconvolute(spectrum, Params(feat), range).ToList();

            Assert.AreEqual(1, envelopes.Count);
            var env = envelopes[0];
            Assert.AreEqual(feat.Charge, env.Charge);
            Assert.AreEqual(feat.Mz.ToMass(feat.Charge), env.MonoisotopicMass, 1e-9);
            Assert.AreEqual(feat.Intensity, env.TotalIntensity);
            Assert.AreEqual(1, env.Peaks.Count, "synthetic envelope carries a single peak");
        }

        [Test]
        public void DirectAlgorithm_FeatureRtWindowOutsideRange_NoEmit()
        {
            // Feature elutes at [10, 15]; range asks for [20, 21].
            var feat = Feature();
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 20.0, maxRt: 21.0);

            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(feat), range).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void DirectAlgorithm_FeatureRtWindowOverlapsRange_Emits()
        {
            // Feature [10, 15] overlaps requested [14, 18] at [14, 15].
            var feat = Feature();
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 14.0, maxRt: 18.0);

            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(feat), range).ToList();
            Assert.AreEqual(1, envelopes.Count);
        }

        [Test]
        public void DirectAlgorithm_FeatureMzOutsideMzRange_NoEmit()
        {
            // Feature m/z = 600; range asks for [700, 705].
            var feat = Feature();
            var range = new MzRtRange(minMZ: 700.0, maxMZ: 705.0, minRt: 12.0, maxRt: 13.0);

            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(feat), range).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void DirectAlgorithm_FeatureChargeOutsideAssumedRange_NoEmit()
        {
            // Feature charge = 50, params accept only [1, 10].
            var feat = Feature(charge: 50);
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 12.0, maxRt: 13.0);
            var parameters = new FromFileDeconvolutionParameters(new[] { feat }, minCharge: 1, maxCharge: 10);

            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), parameters, range).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void DirectAlgorithm_RequiresMzRtRange_ThrowsOnPlainMzRange()
        {
            // FromFileDeconvolutionAlgorithm cannot deconvolute without an RT window
            // because the features carry RT bounds — Deconvoluter raises
            // ArgumentException before the algorithm even runs.
            var feat = Feature();
            var range = new MzRange(599.0, 601.0);

            Assert.Throws<System.ArgumentException>(
                () => Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(feat), range).ToList());
        }

        [Test]
        public void DirectAlgorithm_MultipleFeaturesInsideRange_EmitsAll()
        {
            var featA = new SingleChargeMs1Feature(599.5, Charge: 2, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 1e5);
            var featB = new SingleChargeMs1Feature(600.5, Charge: 3, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 2e5);
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 12.0, maxRt: 13.0);

            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(featA, featB), range).ToList();
            Assert.AreEqual(2, envelopes.Count);
            CollectionAssert.AreEquivalent(new[] { 2, 3 }, envelopes.Select(e => e.Charge));
        }

        [Test]
        public void DirectAlgorithm_NoFeatures_NoEmit()
        {
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 12.0, maxRt: 13.0);
            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(), range).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        // -------- integration tests via MsDataScan.GetIsolatedMassesAndCharges --
        // These exercise the MM-style path: GIMC auto-upgrades the requested range
        // to an MzRtRange (with rtTolerance = 0.1) using the MS2 scan's own RT,
        // then dispatches through Deconvoluter.

        [Test]
        public void Integration_GetIsolatedMassesAndCharges_RtInsideAndMzInIsolation_EmitsOne()
        {
            // MS2 at RT 12.5, isolating m/z 600 ± 1. Feature m/z 600, RT [10, 15].
            // GIMC builds MzRtRange around the MS2's RT (±0.1); feature window covers
            // that range; m/z is inside the isolation window; one envelope returned.
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), Params(Feature())).ToList();
            Assert.AreEqual(1, envelopes.Count);
            Assert.AreEqual(2, envelopes[0].Charge);
            Assert.AreEqual(600.0.ToMass(2), envelopes[0].MonoisotopicMass, 1e-9);
        }

        [Test]
        public void Integration_GetIsolatedMassesAndCharges_MzOutsideIsolation_NoEmit()
        {
            // Feature at m/z 600 is inside the algorithm's 8.5-padded m/z range
            // [591.5, 608.5] but OUTSIDE the strict IsolationRange [599, 601] -- wait,
            // 600 is inside that. Use 597 instead, which is in the padded range but
            // outside strict isolation.
            var feat = new SingleChargeMs1Feature(597.0, 2, 10, 15, 1e5);
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), Params(feat)).ToList();
            Assert.AreEqual(0, envelopes.Count,
                "post-filter on isolationRange.Contains(peak.mz) rejects features that the padded m/z range admitted but the strict isolation excludes");
        }

        [Test]
        public void Integration_GetIsolatedMassesAndCharges_RtOutsideFeatureWindow_NoEmit()
        {
            // MS2 at RT 100; feature window [10, 15].
            var ms2 = MakeMs2(2, rt: 100.0, isolationMz: 600.0, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), Params(Feature())).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void Integration_GetIsolatedMassesAndCharges_NullIsolation_ReturnsEmpty()
        {
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: null, isolationWidth: null);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), Params(Feature())).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void Integration_GetIsolatedMassesAndCharges_TwoFeaturesIntoSameIsolation_EmitsChimera()
        {
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 4.0);
            var featA = new SingleChargeMs1Feature(599.0, Charge: 2, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 1e5);
            var featB = new SingleChargeMs1Feature(601.0, Charge: 3, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 2e5);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), Params(featA, featB)).ToList();
            Assert.AreEqual(2, envelopes.Count);
            CollectionAssert.AreEquivalent(new[] { 2, 3 }, envelopes.Select(e => e.Charge));
        }

        [Test]
        public void Integration_OneFeatureSpansMultipleMs2Scans_EachScanResolvesIt()
        {
            // Caller iterates MS2 scans; each call to GIMC sees the same feature.
            var msScans = new[]
            {
                MakeMs2(2, rt: 11.5, isolationMz: 600.0, isolationWidth: 2.0),
                MakeMs2(3, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0),
                MakeMs2(4, rt: 13.5, isolationMz: 600.0, isolationWidth: 2.0),
            };
            var parameters = Params(Feature());
            int totalEnvelopes = msScans
                .Select(ms2 => ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), parameters).Count())
                .Sum();
            Assert.AreEqual(3, totalEnvelopes);
        }

        // -------- end-to-end fixture tests --------------------------------------
        // Load a real .ms1.feature file (FlashDeconv or TopFD format) through the
        // production Ms1FeatureFile reader, then drive the algorithm via
        // MsDataScan.GetIsolatedMassesAndCharges with a synthetic MS2 scan tuned to
        // capture exactly one charge state of feature row #1.
        //
        // Both fixtures contain a feature with Mass ≈ 10835.9 Da (charges 7-18 for
        // FlashDeconv, 7-17 for TopFD) over RT ≈ 2375-2402 s. Charge 10 m/z ≈ 1084.60.
        // An MS2 with isolationMz 1084.6 ± 1.0 at RT 2390 selects only that charge.

        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature",
            /*expectedIntensityPositive*/ false)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature",
            /*expectedIntensityPositive*/ true)]
        public void EndToEnd_RealMs1FeatureFile_ResolvesExpectedChargeState(string relativeFixturePath, bool expectedIntensityPositive)
        {
            var fixturePath = Path.Combine(TestContext.CurrentContext.TestDirectory, relativeFixturePath);
            var perChargeFeatures = new Ms1FeatureFile(fixturePath).GetMs1Features().ToList();
            Assert.IsTrue(perChargeFeatures.Count >= 25,
                $"fixture {Path.GetFileName(relativeFixturePath)} unexpectedly small: {perChargeFeatures.Count} per-charge entries");

            var parameters = new FromFileDeconvolutionParameters(perChargeFeatures, minCharge: 1, maxCharge: 60);
            var ms2 = MakeMs2(2, rt: 2390.0, isolationMz: 1084.6, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), parameters).ToList();

            Assert.AreEqual(1, envelopes.Count, "exactly one feature-charge in this isolation window");
            var env = envelopes[0];
            Assert.AreEqual(10, env.Charge);
            Assert.AreEqual(10835.9, env.MonoisotopicMass, 0.1);
            double monoMz = env.MonoisotopicMass.ToMz(env.Charge);
            Assert.IsTrue(monoMz > 1083.6 && monoMz < 1085.6, $"emitted m/z {monoMz} outside isolation window");

            // FlashDeconv format has no Apex_intensity column → Intensity flows through
            // as 0. TopFD has it → positive.
            if (expectedIntensityPositive)
                Assert.IsTrue(env.TotalIntensity > 0, "TopFD-derived envelope should carry a positive intensity");
            else
                Assert.AreEqual(0, env.TotalIntensity, "FlashDeconv-derived envelope has no apex column; intensity is 0");
        }

        [Test]
        public void EndToEnd_RealMs1FeatureFile_NoMs2InWindow_NoEmit()
        {
            var fixturePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature");
            var perChargeFeatures = new Ms1FeatureFile(fixturePath).GetMs1Features().ToList();
            var parameters = new FromFileDeconvolutionParameters(perChargeFeatures, minCharge: 1, maxCharge: 60);

            // MS2 at RT 100000 — far outside any feature's window.
            var ms2 = MakeMs2(2, rt: 100000.0, isolationMz: 1084.6, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), parameters).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }
    }
}
