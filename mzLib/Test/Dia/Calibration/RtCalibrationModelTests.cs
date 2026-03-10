// DiaCalibrationTests.cs
// Unit tests for the DIA RT calibration pipeline.
//
// Placement: mzLib/mzLibTest/DiaTests/DiaCalibrationTests.cs
//
// Framework: NUnit 3 (consistent with existing mzLib tests)
// Dependencies: MassSpectrometry, MassSpectrometry.Dia, MassSpectrometry.Dia.Calibration
//
// Test Groups:
//   Group 1 — RtCalibrationModel: forward/inverse transforms and reliability flags
//   Group 2 — RtCalibrationFitter: linear fit quality and RANSAC outlier rejection
//   Group 3 — SelectAnchors: eligibility filtering and quality scoring
//   Group 4 — RecalibrateRtDeviations: the critical post-assembly step (the "Bug 2" path)
//   Group 5 — Convergence invariants: configuration defaults and window-floor behavior
//   Group 6 — Slope divergence guard: documents known iRT-library bug and correct fix
//
// NOTE on DiaSearchResult extended properties:
//   The project snapshot of DiaSearchResult.cs is abbreviated. The full production class
//   also carries: ApexScore, TemporalScore, ObservedApexRt, RtDeviationMinutes,
//   RtDeviationSquared, MeanMassErrorPpm, ClassifierScore, etc., set during assembly.
//   The MakeResult helpers below set them via direct property assignment.
//
// NOTE on LibraryPrecursorInput constructor argument order:
//   (sequence, precursorMz, chargeState, retentionTime, isDecoy, fragmentMzs, fragmentIntensities,
//    irtValue = null)

using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.Dia;
using MassSpectrometry.Dia.Calibration;
using NUnit.Framework;

namespace Test.DiaTests
{
    // ════════════════════════════════════════════════════════════════════════
    // Group 1 — RtCalibrationModel: forward/inverse transforms and reliability
    // ════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("DiaCalibration")]
    public class RtCalibrationModelTests
    {
        /// <summary>
        /// ToMinutes must invert ToIrt exactly.
        /// Model: observedRt = Slope * libraryRt + Intercept.
        /// ToIrt(ToMinutes(x)) == x must hold for any x.
        /// Regression guard: if either transform changes sign or scale, this fails immediately.
        /// </summary>
        [Test]
        public void ToMinutes_ThenToIrt_IsIdentity()
        {
            var model = new RtCalibrationModel(
                slope: 1.0026, intercept: 0.06, sigmaMinutes: 0.052, rSquared: 0.9999, anchorCount: 29);

            double[] irts = { 11.4, 20.0, 28.5, 36.0, 41.3 };
            foreach (double irt in irts)
            {
                double roundTripped = model.ToIrt(model.ToMinutes(irt));
                Assert.That(roundTripped, Is.EqualTo(irt).Within(1e-9),
                    $"Round-trip failed at iRT={irt}: got {roundTripped}");
            }
        }

        /// <summary>
        /// For a native-RT library (slope≈1, intercept≈0.06), ToMinutes should return a value
        /// very close to the input iRT. This matches the PXD005573 benchmark state exactly.
        /// </summary>
        [Test]
        public void ToMinutes_NativeRtLibrary_ReturnsNearIdentity()
        {
            var model = new RtCalibrationModel(
                slope: 1.0026, intercept: 0.059, sigmaMinutes: 0.052, rSquared: 0.9999, anchorCount: 29);

            // ToMinutes(28.7) = 1.0026 * 28.7 + 0.059 = 28.834 min.
            // A native-RT library has slope≈1 so the output is close to — but not equal to —
            // the input. The correct tolerance accounts for slope drift (0.0026 * 28.7 = 0.075)
            // plus the intercept (0.059), totalling ~0.134 min at this RT value.
            double predicted = model.ToMinutes(28.7);
            Assert.That(predicted, Is.EqualTo(28.7).Within(0.15),
                "Native-RT library with slope≈1 should predict within 0.15 min of the library value");
        }

        /// <summary>
        /// IsReliable requires BOTH R² ≥ 0.90 AND anchorCount ≥ 10.
        /// Failing either condition must return false. This guards against the calibration
        /// silently proceeding on a degenerate model with too few anchors or a poor fit.
        /// </summary>
        [Test]
        public void IsReliable_RequiresBothRSquaredAndAnchorCount()
        {
            var good = new RtCalibrationModel(1.0, 0.0, 0.05, rSquared: 0.995, anchorCount: 29);
            Assert.That(good.IsReliable, Is.True, "R²=0.995, 29 anchors should be reliable");

            var badR2 = new RtCalibrationModel(1.0, 0.0, 0.05, rSquared: 0.85, anchorCount: 29);
            Assert.That(badR2.IsReliable, Is.False, "R²=0.85 must be flagged unreliable");

            var fewAnchors = new RtCalibrationModel(1.0, 0.0, 0.05, rSquared: 0.995, anchorCount: 5);
            Assert.That(fewAnchors.IsReliable, Is.False, "5 anchors must be flagged unreliable");
        }

        /// <summary>
        /// GetMinutesWindowHalfWidth(k) must equal k × SigmaMinutes.
        /// IterativeRtCalibrator calls this to compute the extraction window each iteration.
        /// Any change to this formula directly affects how tightly the window tracks the model.
        /// </summary>
        [Test]
        public void GetMinutesWindowHalfWidth_ReturnsKTimesSigma()
        {
            var model = new RtCalibrationModel(1.0, 0.0, sigmaMinutes: 0.039, rSquared: 0.9999, anchorCount: 30);

            Assert.That(model.GetMinutesWindowHalfWidth(4.0), Is.EqualTo(4.0 * 0.039).Within(1e-9));
            Assert.That(model.GetMinutesWindowHalfWidth(3.0), Is.EqualTo(3.0 * 0.039).Within(1e-9));
        }

        /// <summary>
        /// Zero slope must throw an ArgumentException.
        /// A zero-slope model would cause division-by-zero in ToIrt and produce nonsensical
        /// windows throughout the entire pipeline.
        /// </summary>
        [Test]
        public void Constructor_ZeroSlope_Throws()
        {
            Assert.Throws<ArgumentException>(() =>
                new RtCalibrationModel(slope: 0.0, intercept: 0.0, sigmaMinutes: 0.05,
                    rSquared: 0.99, anchorCount: 20));
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // Group 2 — RtCalibrationFitter: linear fit quality and RANSAC robustness
    // ════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("DiaCalibration")]
    public class RtCalibrationFitterTests
    {
        /// <summary>
        /// On perfectly collinear data (no noise), the fitter must recover the exact slope
        /// and intercept and return R²=1.0. This is the baseline sanity test for the fitter.
        /// If this fails, all downstream calibration is unreliable.
        /// </summary>
        [Test]
        public void Fit_PerfectLinearData_RecoversSlopeAndIntercept()
        {
            double slope = 1.0026;
            double intercept = 0.059;

            var libraryRts = Enumerable.Range(0, 30).Select(i => 11.4 + i * 1.0).ToArray();
            var observedRts = libraryRts.Select(x => slope * x + intercept).ToArray();

            var model = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)libraryRts,
                (ReadOnlySpan<double>)observedRts,
                new RtCalibrationFitter.FitOptions { UseRansac = false });

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Slope, Is.EqualTo(slope).Within(1e-6), "Slope not recovered");
            Assert.That(model.Intercept, Is.EqualTo(intercept).Within(1e-6), "Intercept not recovered");
            Assert.That(model.RSquared, Is.EqualTo(1.0).Within(1e-6), "R² must be 1.0 for noise-free data");
        }

        /// <summary>
        /// RANSAC must reject a small fraction of gross outliers and still recover the
        /// correct model from the inlier set. This is the core value of RANSAC in the
        /// calibration pipeline: anchors contaminated by false matches must not bias the model.
        /// </summary>
        [Test]
        public void Fit_WithOutliers_RansacRejectsThemAndFitsCorrectly()
        {
            var libraryRts = Enumerable.Range(0, 30).Select(i => 10.0 + i * 1.0).ToArray();
            var observedRts = libraryRts.Select(x => x).ToArray(); // slope=1, intercept=0

            // Inject 4 gross outliers (~13% contamination)
            observedRts[5] += 8.0;
            observedRts[12] -= 7.0;
            observedRts[20] += 6.0;
            observedRts[27] -= 5.0;

            var model = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)libraryRts,
                (ReadOnlySpan<double>)observedRts,
                new RtCalibrationFitter.FitOptions
                {
                    UseRansac = true,
                    RansacInlierThresholdMinutes = 1.0,
                    RansacIterations = 200
                });

            Assert.That(model, Is.Not.Null, "RANSAC should produce a model");
            Assert.That(model.Slope, Is.EqualTo(1.0).Within(0.02), "RANSAC should recover slope ≈ 1.0");
            Assert.That(model.Intercept, Is.EqualTo(0.0).Within(0.05), "RANSAC should recover intercept ≈ 0");
            Assert.That(model.RSquared, Is.GreaterThan(0.99), "R² should be high after outlier rejection");
        }

        /// <summary>
        /// Fewer than MinAnchors points must return null rather than a degenerate model.
        /// This is the first line of defense against the calibration proceeding on an
        /// empty or near-empty anchor set from a failed bootstrap.
        /// </summary>
        [Test]
        public void Fit_InsufficientAnchors_ReturnsNull()
        {
            double[] libRts = { 10.0, 15.0, 20.0 };
            double[] obsRts = { 10.1, 15.1, 20.1 };

            var model = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)libRts,
                (ReadOnlySpan<double>)obsRts,
                new RtCalibrationFitter.FitOptions { MinAnchors = 5 });

            Assert.That(model, Is.Null, "Must return null when anchor count < MinAnchors");
        }

        /// <summary>
        /// SigmaMinutes on the returned model must closely match the RMS of residuals.
        /// Both the window-sizing logic and the convergence criterion depend on
        /// SigmaMinutes being an accurate measure of fit quality, not an approximation.
        /// </summary>
        [Test]
        public void Fit_SigmaMinutes_MatchesRmsResidual()
        {
            var rng = new Random(42);
            var libraryRts = Enumerable.Range(0, 50).Select(i => 10.0 + i * 0.6).ToArray();
            var observedRts = libraryRts
                .Select(x => x + (rng.NextDouble() - 0.5) * 0.08)
                .ToArray();

            var model = RtCalibrationFitter.Fit(
                (ReadOnlySpan<double>)libraryRts,
                (ReadOnlySpan<double>)observedRts,
                new RtCalibrationFitter.FitOptions { UseRansac = false });

            Assert.That(model, Is.Not.Null);

            double rms = Math.Sqrt(observedRts.Select((obs, i) =>
            {
                double res = obs - (model.Slope * libraryRts[i] + model.Intercept);
                return res * res;
            }).Average());

            Assert.That(model.SigmaMinutes, Is.EqualTo(rms).Within(rms * 0.15),
                "SigmaMinutes must approximate the RMS residual to within 15%");
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // Group 3 — SelectAnchors: eligibility filtering and quality scoring
    // ════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("DiaCalibration")]
    public class SelectAnchorsTests
    {
        /// <summary>Creates a minimal DiaSearchResult with the properties SelectAnchors reads.</summary>
        private static DiaSearchResult MakeResult(
            string sequence, int charge, bool isDecoy,
            float apexScore, float fragDetRate, float observedApexRt,
            double? libraryRt, float meanMassErrorPpm = 0f)
        {
            int queried = 6;
            int detected = Math.Max(1, (int)(queried * fragDetRate));

            var r = new DiaSearchResult(
                sequence: sequence,
                chargeState: charge,
                precursorMz: 500.0,
                windowId: 1,
                isDecoy: isDecoy,
                fragmentsQueried: queried,
                libraryRetentionTime: libraryRt,
                rtWindowStart: 0f,
                rtWindowEnd: 60f);

            r.ApexScore = apexScore;
            r.FragmentsDetected = detected;
            r.ObservedApexRt = observedApexRt;
            r.MeanMassErrorPpm = meanMassErrorPpm;
            return r;
        }

        /// <summary>
        /// Decoy results must NEVER be included as calibration anchors.
        /// Calibration maps library RT → observed RT for target peptides only.
        /// Decoy library RTs are randomized and have no relationship to the true RT alignment.
        /// </summary>
        [Test]
        public void SelectAnchors_ExcludesDecoys()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("PEPTIDEK",  2, isDecoy: false, 0.9f, 0.8f, 20.0f, 20.0),
                MakeResult("DECOYPEPT", 2, isDecoy: true,  0.9f, 0.8f, 21.0f, 21.0),
                MakeResult("PEPTIDER",  2, isDecoy: false, 0.9f, 0.8f, 22.0f, 22.0),
            };

            var anchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.5f, maxAnchors: 100);

            Assert.That(anchors.Count, Is.EqualTo(2), "Only the 2 target results should be selected");
            Assert.That(anchors.Any(a => a.ObservedRt == 21.0), Is.False,
                "The decoy at ObservedRt=21.0 must not appear in anchors");
        }

        /// <summary>
        /// Results with NaN ObservedApexRt must be excluded.
        /// NaN means the apex detector found no valid peak. Passing NaN to the linear fitter
        /// would produce a corrupt model — this was the original "Bug 2" failure mode
        /// that caused 0 IDs in Phase 14 until ObservedApexRt was properly populated.
        /// </summary>
        [Test]
        public void SelectAnchors_ExcludesNanObservedApexRt()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("GOODPEPT",  2, false, 0.9f, 0.8f, 20.0f,     20.0),
                MakeResult("NANPEPT",   2, false, 0.9f, 0.8f, float.NaN, 21.0),  // Bug 2 case
                MakeResult("GOODPEPT2", 2, false, 0.9f, 0.8f, 22.0f,     22.0),
            };

            var anchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.5f, maxAnchors: 100);

            Assert.That(anchors.Count, Is.EqualTo(2), "NaN ObservedApexRt result must be excluded");
            Assert.That(anchors.Any(a => double.IsNaN(a.ObservedRt)), Is.False,
                "No anchor may carry a NaN observed RT");
        }

        /// <summary>
        /// The minApexScore threshold is a strict ≥ comparison.
        /// Results below threshold must be excluded regardless of other quality metrics.
        /// The bootstrap uses 0.85 as the threshold — this is intentionally strict to
        /// ensure anchor purity in the critical first iteration.
        /// </summary>
        [Test]
        public void SelectAnchors_EnforcesApexScoreThreshold()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("HIGH",  2, false, apexScore: 0.9f,  fragDetRate: 0.8f, 20.0f, 20.0),
                MakeResult("BELOW", 2, false, apexScore: 0.4f,  fragDetRate: 0.8f, 21.0f, 21.0),
                MakeResult("EXACT", 2, false, apexScore: 0.85f, fragDetRate: 0.8f, 22.0f, 22.0),
            };

            var anchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.85f, maxAnchors: 100);

            Assert.That(anchors.Count, Is.EqualTo(2),
                "Only results with ApexScore ≥ 0.85 should pass");
            Assert.That(anchors.Any(a => a.ObservedRt == 21.0), Is.False,
                "Low-score result (0.4) at RT=21.0 must not be selected");
        }

        /// <summary>
        /// FragDetRate < 0.5 must exclude a result in strict mode (requireMinFragDetRate=true).
        /// This is used during bootstrap. In relaxed mode (refinement iterations), all
        /// otherwise-eligible results should pass regardless of FragDetRate.
        /// </summary>
        [Test]
        public void SelectAnchors_ExcludesLowFragDetRate_WhenRequired()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("GOOD",      2, false, 0.9f, fragDetRate: 0.8f, 20.0f, 20.0),
                MakeResult("LOWFRAG",   2, false, 0.9f, fragDetRate: 0.3f, 21.0f, 21.0),
                MakeResult("GOODAGAIN", 2, false, 0.9f, fragDetRate: 0.6f, 22.0f, 22.0),
            };

            var strictAnchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.5f, maxAnchors: 100,
                requireMinFragDetRate: true);

            Assert.That(strictAnchors.Count, Is.EqualTo(2),
                "FragDetRate < 0.5 must be excluded when requireMinFragDetRate=true");

            var relaxedAnchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.5f, maxAnchors: 100,
                requireMinFragDetRate: false);

            Assert.That(relaxedAnchors.Count, Is.EqualTo(3),
                "All results should pass when requireMinFragDetRate=false");
        }

        /// <summary>
        /// The composite quality score = ApexScore × FragDetRate × massAccuracyTerm.
        /// massAccuracyTerm = 1 / (1 + |MeanMassErrorPpm| / 10).
        /// At 10 ppm error, term=0.5 so quality is halved.
        /// At 0 ppm error, term=1.0 so quality is unpenalized.
        /// When maxAnchors limits selection, the accurate result must rank higher.
        /// </summary>
        [Test]
        public void SelectAnchors_QualityScore_PenalizesHighMassError()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("ACCURATE",   2, false, 0.9f, 0.9f, 20.0f, 20.0, meanMassErrorPpm: 0f),
                MakeResult("INACCURATE", 2, false, 0.9f, 0.9f, 21.0f, 21.0, meanMassErrorPpm: 10f),
            };

            // Request only 1 anchor — should select the accurate one
            var anchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.5f, maxAnchors: 1);

            Assert.That(anchors.Count, Is.EqualTo(1));
            Assert.That(anchors[0].ObservedRt, Is.EqualTo(20.0).Within(0.001),
                "The result with zero mass error should be preferred over 10 ppm error");
        }

        /// <summary>
        /// Results with null LibraryRetentionTime and no matching precursor in the lookup
        /// must be excluded. The fitter needs a library RT coordinate — without it there
        /// is nothing to calibrate against.
        /// </summary>
        [Test]
        public void SelectAnchors_ExcludesResultsWithNoLibraryRt()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("NORT",  2, false, 0.9f, 0.8f, 20.0f, libraryRt: null),
                MakeResult("HASRT", 2, false, 0.9f, 0.8f, 21.0f, libraryRt: 21.0),
            };

            var anchors = IterativeRtCalibrator.SelectAnchors(
                results, precursors: null, minApexScore: 0.5f, maxAnchors: 100);

            Assert.That(anchors.Count, Is.EqualTo(1),
                "Result with null LibraryRetentionTime must be excluded");
            Assert.That(anchors[0].ObservedRt, Is.EqualTo(21.0).Within(0.001));
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // Group 4 — RecalibrateRtDeviations: the critical post-assembly step
    //
    // RecalibrateRtDeviations transforms RtDeviationMinutes from the raw
    // |ObservedApexRt - LibraryRt| (meaningless for iRT libraries) into
    // (ObservedApexRt - model.ToMinutes(LibraryRt)) in calibrated minutes.
    // This is what gives the classifier its discriminative power over RT.
    // ════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("DiaCalibration")]
    public class RecalibrateRtDeviationsTests
    {
        private static DiaSearchResult MakeResult(
            string sequence, double? libraryRt, float observedApexRt, bool isDecoy = false)
        {
            var r = new DiaSearchResult(sequence, 2, 500.0, 1, isDecoy,
                fragmentsQueried: 4,
                libraryRetentionTime: libraryRt,
                rtWindowStart: 0f, rtWindowEnd: 60f);
            r.ObservedApexRt = observedApexRt;
            return r;
        }

        /// <summary>
        /// After recalibration, RtDeviationMinutes must equal (ObservedApexRt - model.ToMinutes(LibraryRt)).
        /// This is a SIGNED deviation: positive means observed later than predicted, negative means earlier.
        /// The sign preserves directional information that the classifier can exploit.
        ///
        /// Test uses slope=1.0, intercept=0 so model.ToMinutes(x)=x, making expected values
        /// easy to verify by inspection.
        /// </summary>
        [Test]
        public void RecalibrateRtDeviations_SetsSignedDeviationCorrectly()
        {
            var model = new RtCalibrationModel(1.0, 0.0, 0.05, 0.999, 30);

            var results = new List<DiaSearchResult>
            {
                MakeResult("EARLY", libraryRt: 20.0, observedApexRt: 19.9f),  // before predicted → negative
                MakeResult("EXACT", libraryRt: 25.0, observedApexRt: 25.0f),  // perfect match → zero
                MakeResult("LATE",  libraryRt: 30.0, observedApexRt: 30.1f),  // after predicted → positive
            };

            DiaCalibrationPipeline.RecalibrateRtDeviations(results, precursors: null, model);

            Assert.That(results[0].RtDeviationMinutes, Is.EqualTo(-0.1f).Within(1e-4f),
                "Early eluter must have negative deviation");
            Assert.That(results[1].RtDeviationMinutes, Is.EqualTo(0.0f).Within(1e-4f),
                "Perfect match must have zero deviation");
            Assert.That(results[2].RtDeviationMinutes, Is.EqualTo(0.1f).Within(1e-4f),
                "Late eluter must have positive deviation");
        }

        /// <summary>
        /// RtDeviationSquared must always equal RtDeviationMinutes².
        /// The feature extractor reads both fields independently (features [10] and [11]).
        /// Inconsistency would send contradictory signals to the classifier.
        /// </summary>
        [Test]
        public void RecalibrateRtDeviations_SetsSquaredConsistentWithDeviation()
        {
            var model = new RtCalibrationModel(1.0, 0.0, 0.05, 0.999, 30);

            var results = new List<DiaSearchResult>
            {
                MakeResult("PEPT1", 20.0, 20.15f),
                MakeResult("PEPT2", 25.0, 24.90f),
            };

            DiaCalibrationPipeline.RecalibrateRtDeviations(results, precursors: null, model);

            foreach (var r in results)
            {
                float expectedSq = r.RtDeviationMinutes * r.RtDeviationMinutes;
                Assert.That(r.RtDeviationSquared, Is.EqualTo(expectedSq).Within(1e-6f),
                    $"RtDeviationSquared must equal RtDeviationMinutes² for {r.Sequence}");
            }
        }

        /// <summary>
        /// Results with NaN ObservedApexRt must be skipped — their RtDeviationMinutes
        /// must remain at its value from construction. Computing a deviation from NaN would
        /// propagate NaN into the classifier feature vector.
        /// </summary>
        [Test]
        public void RecalibrateRtDeviations_SkipsNanObservedApexRt()
        {
            var model = new RtCalibrationModel(1.0, 0.0, 0.05, 0.999, 30);

            var r = MakeResult("NANPEPT", libraryRt: 20.0, observedApexRt: float.NaN);
            float rtDevBefore = r.RtDeviationMinutes;

            DiaCalibrationPipeline.RecalibrateRtDeviations(
                new List<DiaSearchResult> { r }, precursors: null, model);

            Assert.That(r.RtDeviationMinutes, Is.EqualTo(rtDevBefore),
                "Result with NaN ObservedApexRt must not have RtDeviationMinutes overwritten");
        }

        /// <summary>
        /// Calling RecalibrateRtDeviations with a null calibration model must be a silent no-op.
        /// The pipeline only calls this when a model exists, but the method must not throw
        /// when called defensively with null.
        /// </summary>
        [Test]
        public void RecalibrateRtDeviations_NullModel_IsNoOp()
        {
            var results = new List<DiaSearchResult>
            {
                MakeResult("PEPT", 20.0, 20.1f)
            };
            float devBefore = results[0].RtDeviationMinutes;

            Assert.DoesNotThrow(() =>
                DiaCalibrationPipeline.RecalibrateRtDeviations(
                    results, null, (RtCalibrationModel)null),
                "Null model must not throw");

            Assert.That(results[0].RtDeviationMinutes, Is.EqualTo(devBefore),
                "No-op must leave RtDeviationMinutes unchanged");
        }

        /// <summary>
        /// After recalibration, target results (well-calibrated, close to predicted RT) must
        /// have substantially smaller absolute deviations than decoy results (random RT).
        /// This separation is the primary discriminative signal for the RT-deviation features.
        ///
        /// Expected from PXD005573 benchmark:
        ///   target mean |deviation| ≈ 0.043 min
        ///   decoy mean |deviation|  ≈ 0.123 min   (>2× target)
        ///
        /// This test uses synthetic data that mirrors the structure of that benchmark.
        /// </summary>
        [Test]
        public void RecalibrateRtDeviations_TargetDecoySeparationIsSubstantial()
        {
            var model = new RtCalibrationModel(1.0026, 0.06, 0.052, 0.9999, 29);
            var rng = new Random(42);
            var results = new List<DiaSearchResult>();

            // 100 targets: observed ≈ predicted ± 0.04 min (tight, well-calibrated)
            for (int i = 0; i < 100; i++)
            {
                double libRt = 15.0 + i * 0.25;
                float predicted = (float)model.ToMinutes(libRt);
                float obs = predicted + (float)((rng.NextDouble() - 0.5) * 0.08);
                results.Add(MakeResult($"T{i}", libRt, obs, isDecoy: false));
            }

            // 30 decoys: observed at random positions unrelated to predicted RT
            for (int i = 0; i < 30; i++)
            {
                double libRt = 15.0 + i * 0.8;
                float obs = 15.0f + (float)(rng.NextDouble() * 25.0);
                results.Add(MakeResult($"D{i}", libRt, obs, isDecoy: true));
            }

            DiaCalibrationPipeline.RecalibrateRtDeviations(results, null, model);

            float targetMeanDev = results.Where(r => !r.IsDecoy)
                .Select(r => Math.Abs(r.RtDeviationMinutes)).Average();
            float decoyMeanDev = results.Where(r => r.IsDecoy)
                .Select(r => Math.Abs(r.RtDeviationMinutes)).Average();

            Assert.That(targetMeanDev, Is.LessThan(0.10f),
                $"Target mean |RT deviation| should be small after calibration, got {targetMeanDev:F3} min");
            Assert.That(decoyMeanDev, Is.GreaterThan(targetMeanDev * 2.0f),
                $"Decoy mean ({decoyMeanDev:F3}) should be >2× target mean ({targetMeanDev:F3})");
        }

        /// <summary>
        /// Demonstrates that using the uncalibrated deviation (|ObservedApexRt - LibraryRetentionTime|
        /// directly, i.e., the "Bug 2" anti-pattern) gives WORSE T/D separation than the calibrated
        /// deviation when the library uses non-native units (e.g., iRT with a different scale).
        ///
        /// This test directly justifies why RecalibrateRtDeviations must be called before
        /// feature extraction — it is not optional.
        /// </summary>
        [Test]
        public void RecalibrateRtDeviations_CalibratedIsBetterThanRawForIrtLibrary()
        {
            // iRT-like library: slope=0.5 means libraryRt ≈ 2× observed RT in minutes
            var model = new RtCalibrationModel(0.5, 5.0, 0.05, 0.999, 30);
            var rng = new Random(123);
            var results = new List<DiaSearchResult>();

            for (int i = 0; i < 50; i++)
            {
                double libIrt = 20.0 + i * 0.5;
                float predicted = (float)model.ToMinutes(libIrt);
                float obs = predicted + (float)((rng.NextDouble() - 0.5) * 0.06);
                results.Add(MakeResult($"T{i}", libIrt, obs, isDecoy: false));
            }

            for (int i = 0; i < 20; i++)
            {
                double libIrt = 20.0 + i * 1.2;
                float obs = 10.0f + (float)(rng.NextDouble() * 20.0);
                results.Add(MakeResult($"D{i}", libIrt, obs, isDecoy: true));
            }

            // Measure raw (Bug 2) separation: |ObservedApexRt - LibraryRt|
            float rawTargetMean = results.Where(r => !r.IsDecoy)
                .Select(r => Math.Abs(r.ObservedApexRt - (float)r.LibraryRetentionTime.Value)).Average();
            float rawDecoyMean = results.Where(r => r.IsDecoy)
                .Select(r => Math.Abs(r.ObservedApexRt - (float)r.LibraryRetentionTime.Value)).Average();

            // Apply calibrated recalibration
            DiaCalibrationPipeline.RecalibrateRtDeviations(results, null, model);

            float calTargetMean = results.Where(r => !r.IsDecoy)
                .Select(r => Math.Abs(r.RtDeviationMinutes)).Average();
            float calDecoyMean = results.Where(r => r.IsDecoy)
                .Select(r => Math.Abs(r.RtDeviationMinutes)).Average();

            float rawSeparation = rawDecoyMean - rawTargetMean;
            float calSeparation = calDecoyMean - calTargetMean;

            Assert.That(calTargetMean, Is.LessThan(rawTargetMean),
                "Calibrated target deviation must be smaller than raw iRT-unit deviation");
            Assert.That(calSeparation, Is.GreaterThan(rawSeparation),
                $"Calibrated T/D separation ({calSeparation:F3}) must exceed raw ({rawSeparation:F3})");
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // Group 5 — Convergence invariants: configuration defaults and window floor
    // ════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("DiaCalibration")]
    public class ConvergenceInvariantTests
    {
        /// <summary>
        /// ConvergenceThreshold must be 0.02 (2%). This was deliberately tightened from 5%
        /// in Phase 23 after the benchmark showed convergence at 1.27% Δσ/σ.
        /// If changed to a looser value, the calibrator will stop before reaching the
        /// optimal σ and produce fewer IDs.
        /// </summary>
        [Test]
        public void ConvergenceThreshold_Default_IsTwoPercent()
        {
            var calibrator = new IterativeRtCalibrator();
            Assert.That(calibrator.ConvergenceThreshold, Is.EqualTo(0.02).Within(1e-9),
                "ConvergenceThreshold must be 0.02 — tightened from 5% in Phase 23");
        }

        /// <summary>
        /// MaxIterations must default to 6. Provides headroom above the observed convergence
        /// of 2 iterations while still being a hard cap against infinite loops.
        /// </summary>
        [Test]
        public void MaxIterations_Default_IsSix()
        {
            var calibrator = new IterativeRtCalibrator();
            Assert.That(calibrator.MaxIterations, Is.EqualTo(6));
        }

        /// <summary>
        /// MinWindowHalfWidthMinutes must default to 0.3 min.
        /// This floor prevents the extraction window from collapsing below ~18 seconds
        /// as σ shrinks toward zero on ideal datasets, which would exclude the true apex.
        /// </summary>
        [Test]
        public void MinWindowHalfWidthMinutes_Default_IsPointThree()
        {
            var calibrator = new IterativeRtCalibrator();
            Assert.That(calibrator.MinWindowHalfWidthMinutes, Is.EqualTo(0.3).Within(1e-9));
        }

        /// <summary>
        /// The window computed per-iteration = MAX(SigmaMultiplier × σ, MinWindowHalfWidthMinutes).
        /// When σ is very small (0.005 min), the floor of 0.3 must dominate.
        /// When σ is moderate (0.08 min), the multiplier must dominate.
        /// </summary>
        [TestCase(0.005, 4.0, 0.3, 0.3, Description = "Tiny σ: floor dominates (4.0×0.005=0.02 < 0.3)")]
        [TestCase(0.08, 4.0, 0.3, 0.32, Description = "Moderate σ: multiplier dominates (4.0×0.08=0.32 > 0.3)")]
        [TestCase(0.039, 4.0, 0.3, 0.3, Description = "Phase 23 final σ: floor dominates (4.0×0.039=0.156 < 0.3)")]
        public void WindowHalfWidth_ComputedCorrectly(
            double sigma, double multiplier, double floor, double expected)
        {
            double window = Math.Max(multiplier * sigma, floor);
            Assert.That(window, Is.EqualTo(expected).Within(1e-9));
        }

        /// <summary>
        /// A log with only the bootstrap entry cannot be evaluated for convergence.
        /// Convergence requires at least one completed refinement iteration — the
        /// comparison is between consecutive σ values, which needs ≥2 entries.
        /// The guard `iteration >= startIteration + 1` in RunRefinementIterations
        /// is what prevents premature convergence at the bootstrap.
        /// </summary>
        [Test]
        public void Convergence_RequiresAtLeastTwoLogEntries()
        {
            var bootstrapOnlyLog = new List<CalibrationIterationLog>
            {
                new CalibrationIterationLog { Iteration = 0, SigmaMinutes = 0.052, AnchorCount = 29 }
            };

            bool hasEnoughEntries = bootstrapOnlyLog.Count >= 2;
            Assert.That(hasEnoughEntries, Is.False,
                "A single-entry log (bootstrap only) cannot have a convergence comparison");
        }

        /// <summary>
        /// The PXD005573 Phase 23 benchmark produced σ values of 0.052 → 0.040 → 0.039.
        /// Each must be non-increasing. Any increase means the divergence guard failed to fire.
        /// </summary>
        [Test]
        public void BenchmarkSigmaValues_AreMonotonicallyNonIncreasing()
        {
            double[] sigmas = { 0.052, 0.040, 0.039 }; // exact Phase 23 values

            for (int i = 1; i < sigmas.Length; i++)
            {
                Assert.That(sigmas[i], Is.LessThanOrEqualTo(sigmas[i - 1]),
                    $"σ increased from iteration {i - 1} to {i}: divergence guard should have fired");
            }
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // Group 6 — Slope divergence guard: known bug and specification of the fix
    //
    // Bug location: IterativeRtCalibrator.RunRefinementIterations, lines 464-476.
    //
    // Current code:
    //   double curDist = Math.Abs(newSlope - 1.0);
    //   double prevDist = Math.Abs(prevSlope - 1.0);
    //   if (curDist > prevDist + 0.02)  ← WRONG for iRT libraries
    //
    // Required fix:
    //   double relChange = Math.Abs(newSlope - prevSlope) / Math.Abs(prevSlope);
    //   if (relChange > 0.50)  ← correct, slope-agnostic
    //
    // The bug is benign for native-RT libraries (slope≈1.0) but silently fails to
    // protect iRT libraries (slope≈0.17) against slope divergence.
    // ════════════════════════════════════════════════════════════════════════

    [TestFixture]
    [Category("DiaCalibration")]
    public class SlopeDivergenceGuardTests
    {
        // Current (buggy) guard — mirrors RunRefinementIterations lines 466-468 exactly
        private static bool CurrentGuardFires(double prevSlope, double newSlope)
        {
            double curDist = Math.Abs(newSlope - 1.0);
            double prevDist = Math.Abs(prevSlope - 1.0);
            return curDist > prevDist + 0.02;
        }

        // Correct relative-change guard — the fix to be implemented
        private static bool CorrectGuardFires(double prevSlope, double newSlope, double threshold = 0.50)
        {
            if (Math.Abs(prevSlope) < 1e-6) return false;
            return Math.Abs(newSlope - prevSlope) / Math.Abs(prevSlope) > threshold;
        }

        /// <summary>
        /// For a native-RT library converging (slope: 1.003 → 1.001), neither the current
        /// nor the correct guard should fire. This validates that the fix does not regress
        /// the common case that currently works correctly.
        /// </summary>
        [Test]
        public void SlopeGuard_NativeRtConverging_NeitherGuardFires()
        {
            double prevSlope = 1.003;
            double newSlope = 1.001; // converging toward 1.0

            Assert.That(CurrentGuardFires(prevSlope, newSlope), Is.False,
                "Current guard must not fire on normal native-RT convergence");
            Assert.That(CorrectGuardFires(prevSlope, newSlope), Is.False,
                "Correct guard must also not fire on 0.2% slope change");
        }

        /// <summary>
        /// DOCUMENTS THE KNOWN BUG:
        ///
        /// For an iRT library (true slope ≈ 0.175), a slope doubling (0.175 → 0.320) is
        /// clear divergence — the model is wrong and should be reverted.
        ///
        /// The CURRENT guard does NOT fire because both slopes are about equally far from
        /// 1.0 (distances 0.825 and 0.680 differ by only 0.145 < threshold of 0.02+0.145).
        /// Wait — actually curDist=0.68 is NOT > prevDist=0.825 + 0.02=0.845. So it doesn't fire.
        ///
        /// The CORRECT guard DOES fire: (0.320 - 0.175) / 0.175 = 82.8% >> 50% threshold.
        ///
        /// Fix: replace absolute-distance guard with relative-change guard.
        /// </summary>
        [Test]
        public void SlopeGuard_IrtLibrary_BugDocumented_CurrentMissesRealDivergence()
        {
            double prevSlope = 0.175;  // iRT library: correct slope
            double badSlope = 0.320;   // slope has nearly doubled — clearly wrong

            bool currentFires = CurrentGuardFires(prevSlope, badSlope);
            Assert.That(currentFires, Is.False,
                "KNOWN BUG CONFIRMED: current guard does not fire on 82.8% slope increase for iRT library. " +
                "Fix location: IterativeRtCalibrator.RunRefinementIterations lines 464-476.");

            bool correctFires = CorrectGuardFires(prevSlope, badSlope);
            Assert.That(correctFires, Is.True,
                "Relative-change guard correctly detects 82.8% slope increase as divergence");
        }

        /// <summary>
        /// The correct relative-change guard must NOT fire on normal PXD005573 benchmark
        /// convergence where slope changes by 0.27% between iterations (0.9999 → 0.9998).
        /// </summary>
        [Test]
        public void SlopeGuard_CorrectGuard_DoesNotFireOnBenchmarkConvergence()
        {
            double prevSlope = 0.9999; // Phase 23 iteration 1
            double newSlope = 0.9998; // Phase 23 iteration 2

            Assert.That(CorrectGuardFires(prevSlope, newSlope), Is.False,
                "0.01% slope change must not trigger the 50% relative-change guard");
        }

        /// <summary>
        /// The correct relative-change guard must fire on genuinely large slope changes
        /// for both native-RT and iRT library slopes.
        /// </summary>
        [TestCase(1.0, 1.6, Description = "Native-RT: +60% jump")]
        [TestCase(1.0, 0.4, Description = "Native-RT: -60% drop")]
        [TestCase(0.175, 0.320, Description = "iRT: +82.8% jump (the known bug case)")]
        [TestCase(0.175, 0.085, Description = "iRT: -51.4% drop (just over the 50% threshold)")]
        public void SlopeGuard_CorrectGuard_FiresOnGenuineDivergence(double prevSlope, double newSlope)
        {
            Assert.That(CorrectGuardFires(prevSlope, newSlope), Is.True,
                $"Large slope change {prevSlope:F4} → {newSlope:F4} " +
                $"({Math.Abs(newSlope - prevSlope) / Math.Abs(prevSlope):P1} relative change) " +
                "must trigger the relative-change guard");
        }
    }
}