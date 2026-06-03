using System;
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

        // Tests use the internal IEnumerable ctor (InternalsVisibleTo("Test") is declared in
        // Readers) to seed parameters without writing a temp feature file.
        private static FromFileDeconvolutionParameters Params(params ISingleChargeMs1Feature[] feats)
            => new FromFileDeconvolutionParameters(feats, minCharge: 1, maxCharge: 60);

        // -------- direct algorithm tests via Deconvoluter.Deconvolute -----------
        // These tests bypass MsDataScan and exercise FromFileDeconvolutionAlgorithm
        // through Deconvoluter directly, with hand-built MzRtRange inputs.

        [Test]
        public void DirectAlgorithm_FeatureInsideRange_EmitsOneEnvelope()
        {
            var feat = Feature();
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 12.0, maxRt: 13.0);

            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(feat), range).ToList();

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
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 20.0, maxRt: 21.0);
            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(Feature()), range).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void DirectAlgorithm_FeatureRtWindowOverlapsRange_Emits()
        {
            // Feature [10, 15] overlaps requested [14, 18] at [14, 15].
            var range = new MzRtRange(minMZ: 599.0, maxMZ: 601.0, minRt: 14.0, maxRt: 18.0);
            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(Feature()), range).ToList();
            Assert.AreEqual(1, envelopes.Count);
        }

        [Test]
        public void DirectAlgorithm_FeatureMzOutsideMzRange_NoEmit()
        {
            // Feature m/z = 600; range asks for [700, 705].
            var range = new MzRtRange(minMZ: 700.0, maxMZ: 705.0, minRt: 12.0, maxRt: 13.0);
            var envelopes = Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(Feature()), range).ToList();
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
            // FromFile requires an MzRtRange. Deconvoluter raises ArgumentException
            // before the algorithm even runs.
            var range = new MzRange(599.0, 601.0);
            Assert.Throws<System.ArgumentException>(
                () => Deconvoluter.Deconvolute(MinimalMs1Spectrum(), Params(Feature()), range).ToList());
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

        // -------- contract tests for the MzSpectrum / MsDataScan overloads ------

        [Test]
        public void MzSpectrumOverloadOfGimc_FromFileWithoutMzRtRange_Throws()
        {
            // MsDataScan.GetIsolatedMassesAndCharges(MzSpectrum, ...) builds a plain
            // MzRange and hands off to Deconvoluter.Deconvolute(MzSpectrum, ...),
            // which throws when the params demand RT context. Documenting that
            // contract so callers know to use the MsDataScan overload instead.
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0);
            Assert.Throws<System.ArgumentException>(
                () => ms2.GetIsolatedMassesAndCharges(MinimalMs1Spectrum(), Params(Feature())).ToList());
        }

        // -------- integration tests via the MsDataScan overload of GIMC ---------
        // The MM-side flow: pass a precursor MS1 MsDataScan; Deconvoluter auto-upgrades
        // the plain MzRange to an MzRtRange using the precursor's RT.

        [Test]
        public void Integration_GimcMsDataScanOverload_RtInsideAndMzInIsolation_EmitsOne()
        {
            var precursorMs1 = MakeMs1(1, rt: 12.5);
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, Params(Feature())).ToList();
            Assert.AreEqual(1, envelopes.Count);
            Assert.AreEqual(2, envelopes[0].Charge);
            Assert.AreEqual(600.0.ToMass(2), envelopes[0].MonoisotopicMass, 1e-9);
        }

        [Test]
        public void Integration_GimcMsDataScanOverload_MzOutsideIsolation_NoEmit()
        {
            // m/z 597 is inside the 8.5-padded m/z range [591.5, 608.5] but OUTSIDE
            // the strict isolation [599, 601]. The post-filter rejects it.
            var feat = new SingleChargeMs1Feature(597.0, 2, 10, 15, 1e5);
            var precursorMs1 = MakeMs1(1, rt: 12.5);
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, Params(feat)).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void Integration_GimcMsDataScanOverload_RtOutsideFeatureWindow_NoEmit()
        {
            var precursorMs1 = MakeMs1(1, rt: 100.0);
            var ms2 = MakeMs2(2, rt: 100.0, isolationMz: 600.0, isolationWidth: 2.0);
            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, Params(Feature())).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void Integration_GimcMsDataScanOverload_NullIsolation_ReturnsEmpty()
        {
            var precursorMs1 = MakeMs1(1, rt: 12.5);
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: null, isolationWidth: null);
            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, Params(Feature())).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        [Test]
        public void Integration_GimcMsDataScanOverload_TwoFeaturesIntoSameIsolation_EmitsChimera()
        {
            var precursorMs1 = MakeMs1(1, rt: 12.5);
            var ms2 = MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 4.0);
            var featA = new SingleChargeMs1Feature(599.0, Charge: 2, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 1e5);
            var featB = new SingleChargeMs1Feature(601.0, Charge: 3, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 2e5);

            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, Params(featA, featB)).ToList();
            Assert.AreEqual(2, envelopes.Count);
            CollectionAssert.AreEquivalent(new[] { 2, 3 }, envelopes.Select(e => e.Charge));
        }

        [Test]
        public void Integration_OneFeatureSpansMultipleMs2Scans_EachScanResolvesIt()
        {
            var parameters = Params(Feature());
            var msPairs = new (MsDataScan precursor, MsDataScan ms2)[]
            {
                (MakeMs1(1, rt: 11.5), MakeMs2(2, rt: 11.5, isolationMz: 600.0, isolationWidth: 2.0)),
                (MakeMs1(3, rt: 12.5), MakeMs2(4, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0)),
                (MakeMs1(5, rt: 13.5), MakeMs2(6, rt: 13.5, isolationMz: 600.0, isolationWidth: 2.0)),
            };
            int totalEnvelopes = msPairs
                .Select(p => p.ms2.GetIsolatedMassesAndCharges(p.precursor, parameters).Count())
                .Sum();
            Assert.AreEqual(3, totalEnvelopes);
        }

        // -------- end-to-end fixture tests --------------------------------------
        // Drive the file-path constructor for FromFileDeconvolutionParameters,
        // exercising the format-detection path (FileReader.ReadResultFile +
        // SupportedFileType + Ms1FeatureFile.LoadResults). Each fixture loads a real
        // .ms1.feature file (FlashDeconv or TopFD), then pairs against a synthetic
        // MS2 tuned to capture exactly one charge state of feature row #1.

        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature",
            /*expectedIntensityPositive*/ false)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature",
            /*expectedIntensityPositive*/ true)]
        public void EndToEnd_FilePathCtor_ResolvesExpectedChargeState(string relativeFixturePath, bool expectedIntensityPositive)
        {
            var fixturePath = Path.Combine(TestContext.CurrentContext.TestDirectory, relativeFixturePath);

            // Drives the public file-path ctor — exercises FileReader.ReadResultFile,
            // SupportedFileType auto-detection, and Ms1FeatureFile.LoadResults.
            var parameters = new FromFileDeconvolutionParameters(fixturePath, minCharge: 1, maxCharge: 60);
            Assert.IsTrue(parameters.Features.Count >= 25,
                $"fixture {Path.GetFileName(relativeFixturePath)} unexpectedly small: {parameters.Features.Count} per-charge entries");

            var precursorMs1 = MakeMs1(1, rt: 2390.0);
            var ms2 = MakeMs2(2, rt: 2390.0, isolationMz: 1084.6, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, parameters).ToList();

            Assert.AreEqual(1, envelopes.Count, "exactly one feature-charge in this isolation window");
            var env = envelopes[0];
            Assert.AreEqual(10, env.Charge);
            Assert.AreEqual(10835.9, env.MonoisotopicMass, 0.1);
            double monoMz = env.MonoisotopicMass.ToMz(env.Charge);
            Assert.IsTrue(monoMz > 1083.6 && monoMz < 1085.6, $"emitted m/z {monoMz} outside isolation window");

            if (expectedIntensityPositive)
                Assert.IsTrue(env.TotalIntensity > 0, "TopFD-derived envelope should carry a positive intensity");
            else
                Assert.AreEqual(0, env.TotalIntensity, "FlashDeconv-derived envelope has no apex column; intensity is 0");
        }

        [Test]
        public void EndToEnd_FilePathCtor_NoMs2InWindow_NoEmit()
        {
            var fixturePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature");
            var parameters = new FromFileDeconvolutionParameters(fixturePath, minCharge: 1, maxCharge: 60);

            // MS2 at RT 100000 — far outside any feature's window.
            var precursorMs1 = MakeMs1(1, rt: 100000.0);
            var ms2 = MakeMs2(2, rt: 100000.0, isolationMz: 1084.6, isolationWidth: 2.0);

            var envelopes = ms2.GetIsolatedMassesAndCharges(precursorMs1, parameters).ToList();
            Assert.AreEqual(0, envelopes.Count);
        }

        // -------- contract / error-path tests for FromFileDeconvolutionParameters

        [Test]
        public void FromFileDeconvolutionParameters_NullFilePath_Throws()
        {
            Assert.Throws<System.ArgumentNullException>(
                () => new FromFileDeconvolutionParameters((string)null!, minCharge: 1, maxCharge: 60));
        }

        [Test]
        public void FromFileDeconvolutionParameters_NonFeatureFile_ThrowsMzLibException()
        {
            // .psmtsv is a recognized SupportedFileType but its IResultFile (PsmFromTsvFile)
            // doesn't implement IMs1FeatureFile — the ctor must reject it with a clear
            // message rather than silently producing an empty feature set.
            var psmTsv = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");
            Assert.Throws<MzLibUtil.MzLibException>(
                () => new FromFileDeconvolutionParameters(psmTsv, minCharge: 1, maxCharge: 60));
        }

        [Test]
        public void FromFileDeconvolutionParameters_ToDecoyParameters_ReturnsNull()
        {
            // FromFile has no decoy variant (masses are taken as authoritative).
            Assert.IsNull(Params(Feature()).ToDecoyParameters());
        }

        // -------- Deconvoluter auto-upgrade with null range ---------------------

        [Test]
        public void Deconvoluter_MsDataScanOverload_NullRange_AutoUpgradesUsingScanRange()
        {
            // When the caller passes no range, Deconvolute(MsDataScan, ...) for FromFile
            // params auto-upgrades using the scan's MassSpectrum.Range and RetentionTime.
            // Covers the null-side of the ternary in the auto-upgrade branch.
            var precursorMs1 = MakeMs1(1, rt: 12.5);
            // Feature m/z 500 lies inside MinimalMs1Spectrum's range (single peak at 500.0).
            var feat = new SingleChargeMs1Feature(500.0, Charge: 2, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 1e5);

            var envelopes = Deconvoluter.Deconvolute(precursorMs1, Params(feat)).ToList();
            Assert.AreEqual(1, envelopes.Count);
            Assert.AreEqual(2, envelopes[0].Charge);
        }

        [Test]
        public void Deconvoluter_DeconvoluteWithDecoys_FromFileParams_ThrowsBecauseNoDecoySupport()
        {
            // FromFileDeconvolutionParameters.ToDecoyParameters() returns null.
            // DeconvoluteWithDecoys must surface that as InvalidOperationException
            // rather than crashing deeper. Locks Deconvoluter's documented contract.
            var feat = Feature();
            var range = new MzRtRange(599.0, 601.0, 12.0, 13.0);
            Assert.Throws<System.InvalidOperationException>(
                () => Deconvoluter.DeconvoluteWithDecoys(MinimalMs1Spectrum(), Params(feat), range));
        }

        // Test-only stub: a DeconvolutionParameters that claims DeconvolutionType.FromFile
        // but doesn't override CreateAlgorithm. This is the case Deconvoluter's
        // enum-switch FromFile branch is defending against; a future subclass with this
        // shape would otherwise fall through and produce a confusing failure.
        private sealed class FromFileParamsWithoutAlgorithmOverride : DeconvolutionParameters
        {
            public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.FromFile;
            public FromFileParamsWithoutAlgorithmOverride() : base(minCharge: 1, maxCharge: 10) { }
            public override DeconvolutionParameters ToDecoyParameters() => null;
            protected override bool EqualProperties(DeconvolutionParameters other) => true;
            protected override void AddHashCodes(HashCode hash) { }
            public override FromFileParamsWithoutAlgorithmOverride Clone()
            {
                return new FromFileParamsWithoutAlgorithmOverride
                {
                    UseGenericScore = UseGenericScore
                };
            }
        }

        [Test]
        public void Deconvoluter_FromFileWithoutCreateAlgorithmOverride_ThrowsExplanatoryMzLibException()
        {
            var range = new MzRtRange(599.0, 601.0, 12.0, 13.0);
            var brokenParams = new FromFileParamsWithoutAlgorithmOverride();
            var ex = Assert.Throws<MzLibUtil.MzLibException>(
                () => Deconvoluter.Deconvolute(MinimalMs1Spectrum(), brokenParams, range).ToList());
            Assert.IsTrue(ex.Message.Contains("CreateAlgorithm"),
                $"exception message should point at the missing CreateAlgorithm override; was: {ex.Message}");
        }
    }
}
