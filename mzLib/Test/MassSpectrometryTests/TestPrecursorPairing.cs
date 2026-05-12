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
    public sealed class TestPrecursorPairing
    {
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

        [Test]
        public void HappyPath_RtInsideAndMzInIsolation_EmitsOnePair()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var feat = Feature();
            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { feat }, file).ToList();

            Assert.AreEqual(1, pairs.Count);
            var (ms2Scan, env) = pairs[0];
            Assert.AreEqual(feat.Charge, env.Charge);
            Assert.AreEqual(feat.Mz.ToMass(feat.Charge), env.MonoisotopicMass, 1e-9);
            Assert.AreEqual(feat.Intensity, env.TotalIntensity);
            Assert.AreSame(file.GetMsDataScans()[1], ms2Scan);
            Assert.AreEqual(2, ms2Scan.OneBasedScanNumber);
            Assert.AreEqual(1, env.Peaks.Count, "envelope from per-charge feature carries a single synthetic peak");
        }

        [Test]
        public void RtBeforeWindow_NoEmit()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 5.0),
                MakeMs2(2, rt: 9.99, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void RtAfterWindow_NoEmit()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 16.0),
                MakeMs2(2, rt: 15.01, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void MzBelowIsolationWindow_NoEmit()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 700.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature(mz: 600.0) }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void MzAboveIsolationWindow_NoEmit()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 500.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature(mz: 600.0) }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void TwoFeaturesIntoSameMs2_EmitsChimera()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 4.0),
            });

            var featA = new SingleChargeMs1Feature(599.0, Charge: 2, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 1e5);
            var featB = new SingleChargeMs1Feature(601.0, Charge: 3, RetentionTimeStart: 10, RetentionTimeEnd: 15, Intensity: 2e5);

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { featA, featB }, file).ToList();
            Assert.AreEqual(2, pairs.Count);
            Assert.AreEqual(1, pairs.Count(p => p.PrecursorEnvelope.Charge == 2 && p.PrecursorEnvelope.MonoisotopicMass == 599.0.ToMass(2)));
            Assert.AreEqual(1, pairs.Count(p => p.PrecursorEnvelope.Charge == 3 && p.PrecursorEnvelope.MonoisotopicMass == 601.0.ToMass(3)));
        }

        [Test]
        public void OneFeatureSpansMultipleMs2Scans_EmitsOnePerScan()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 11.0),
                MakeMs2(2, rt: 11.5, isolationMz: 600.0, isolationWidth: 2.0),
                MakeMs2(3, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0),
                MakeMs2(4, rt: 13.5, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(3, pairs.Count);
            CollectionAssert.AreEquivalent(new[] { 2, 3, 4 }, pairs.Select(p => p.Ms2Scan.OneBasedScanNumber));
        }

        [Test]
        public void Ms2WithNullIsolation_Skipped()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: null, isolationWidth: null),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void Ms1ScansAreNotPairingTargets()
        {
            // MS1 at RT 12.5 with the feature m/z in its scan range — still must NOT pair
            // because the join is restricted to MS2 scans.
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.5),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void Ms3ScansAreNotPairingTargets()
        {
            // MS3 isolation typically targets an MS2 fragment, not an MS1 precursor.
            // The join restricts pairing to MsnOrder == 2.
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0, msnOrder: 3),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void FeatureRtBoundariesAreInclusive()
        {
            // Scan at exactly RTStart and exactly RTEnd should both match.
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 9.0),
                MakeMs2(2, rt: 10.0, isolationMz: 600.0, isolationWidth: 2.0),
                MakeMs2(3, rt: 15.0, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(2, pairs.Count);
        }

        [Test]
        public void IsolationWindowBoundariesAreInclusive()
        {
            // m/z at exactly Isolation min and exactly Isolation max should both match.
            // Window centered at 600 with width 2.0 -> [599.0, 601.0].
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var featLow = new SingleChargeMs1Feature(599.0, 2, 10, 15, 1e5);
            var featHigh = new SingleChargeMs1Feature(601.0, 2, 10, 15, 1e5);

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { featLow, featHigh }, file).ToList();
            Assert.AreEqual(2, pairs.Count);
        }

        [Test]
        public void NoFeatures_NoEmit()
        {
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 12.0),
                MakeMs2(2, rt: 12.5, isolationMz: 600.0, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new SingleChargeMs1Feature[0], file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        [Test]
        public void NoMs2Scans_NoEmit()
        {
            // A file containing only MS1 scans (no fragmentation) yields no pairs
            // regardless of how many features fall in any RT window.
            var file = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 11.0),
                MakeMs1(2, rt: 12.0),
                MakeMs1(3, rt: 13.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(new[] { Feature() }, file).ToList();
            Assert.AreEqual(0, pairs.Count);
        }

        // End-to-end fixture tests: load a real .ms1.feature file (FlashDeconv or TopFD
        // format) through the real Ms1FeatureFile reader, then pair against a synthetic
        // MS2 scan tuned to capture exactly one charge state of feature row #1.
        //
        // Both fixtures contain a feature with Mass ≈ 10835.9 Da spanning charges 7-18
        // (FlashDeconv) or 7-17 (TopFD) over a similar RT window (≈2375-2402 s).
        // Charge 10 m/z ≈ 1084.60. An isolation window of 1084.6 ± 1.0 catches that
        // single charge: charge 9 m/z ≈ 1204.99, charge 11 m/z ≈ 986.08 — both outside.
        // No other feature in either fixture has a per-charge entry inside [1083.6, 1085.6]
        // at RT 2390, so the join is unambiguous.
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature",
            /*expectedIntensityPositive*/ false)]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature",
            /*expectedIntensityPositive*/ true)]
        public void EndToEnd_RealMs1FeatureFile_PairsExpectedChargeState(string relativeFixturePath, bool expectedIntensityPositive)
        {
            // Load the fixture through the real production reader, exactly as a search
            // pipeline would.
            var fixturePath = Path.Combine(TestContext.CurrentContext.TestDirectory, relativeFixturePath);
            var ms1FeatureFile = new Ms1FeatureFile(fixturePath);
            var perChargeFeatures = ms1FeatureFile.GetMs1Features().ToList();

            // Sanity: the fixtures expand to ≥ 25 per-charge entries (charge ranges plus
            // singly-charged residuals). If this drops, the fixture changed under us.
            Assert.IsTrue(perChargeFeatures.Count >= 25,
                $"fixture {Path.GetFileName(relativeFixturePath)} unexpectedly small: {perChargeFeatures.Count} per-charge entries");

            // Build a single MS2 scan whose isolation window catches only the charge-10
            // expansion of the ≈10835.9 Da feature.
            var ms2File = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 2390.0),
                MakeMs2(2, rt: 2390.0, isolationMz: 1084.6, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(perChargeFeatures, ms2File).ToList();

            Assert.AreEqual(1, pairs.Count, "exactly one feature × MS2 pair expected");
            var (ms2Scan, env) = pairs[0];
            Assert.AreEqual(10, env.Charge);
            Assert.AreEqual(10835.9, env.MonoisotopicMass, 0.1);
            double monoMz = env.MonoisotopicMass.ToMz(env.Charge);
            Assert.IsTrue(monoMz > 1083.6 && monoMz < 1085.6,
                $"emitted m/z {monoMz} is outside the isolation window");
            Assert.AreEqual(2, ms2Scan.OneBasedScanNumber, "must pair to the MS2 scan, not the MS1");

            // FlashDeconv format has no Apex_intensity column → SingleChargeMs1Feature.Intensity
            // is 0. TopFD format does → Intensity is the per-feature reported apex.
            if (expectedIntensityPositive)
                Assert.IsTrue(env.TotalIntensity > 0, "TopFD-derived envelope should carry a positive intensity");
            else
                Assert.AreEqual(0, env.TotalIntensity, "FlashDeconv-derived envelope has no apex column; intensity flows through as 0");
        }

        [Test]
        public void EndToEnd_RealMs1FeatureFile_NoMs2InWindow_NoEmit()
        {
            // Same fixture, but the MS2 scan now sits at an RT far outside any feature's
            // window. The reader and join still run end-to-end; every candidate drops.
            var fixturePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature");
            var perChargeFeatures = new Ms1FeatureFile(fixturePath).GetMs1Features().ToList();

            var ms2File = new FakeMsDataFile(new[]
            {
                MakeMs1(1, rt: 100000.0),
                MakeMs2(2, rt: 100000.0, isolationMz: 1084.6, isolationWidth: 2.0),
            });

            var pairs = Deconvoluter.PairPrecursorsToMs2(perChargeFeatures, ms2File).ToList();
            Assert.AreEqual(0, pairs.Count);
        }
    }
}
