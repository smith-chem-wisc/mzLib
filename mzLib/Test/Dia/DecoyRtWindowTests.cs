// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MassSpectrometry.Dia.Calibration;
using MzLibUtil;
using NUnit.Framework;

namespace Test.DiaTests
{
    /// <summary>
    /// Tests for Option A decoy RT window fix in GenerateFromLowess.
    ///
    /// Root cause of backwards TemporalScore / SpectralAngle / BoundarySignalRatio:
    ///   GenerateFromLowess currently gives decoys the same RT window as their paired
    ///   target, regardless of the decoy's own RetentionTime. Decoys therefore extract
    ///   signal from the target's elution region — real co-eluting peptide signal that
    ///   inflates their chromatographic scores.
    ///
    ///   Current code (fill loop, decoy path):
    ///     int pairedIdx = i - firstDecoyIdx;
    ///     rtMin = targetRtMin[pairedIdx];   // ← ignores decoy's own RT entirely
    ///     rtMax = targetRtMax[pairedIdx];
    ///
    /// Fix:
    ///   Pre-compute decoy windows from each decoy's own RetentionTime through LOWESS,
    ///   exactly as targets do. Fall back to paired target window only when the decoy
    ///   has no RetentionTime.
    ///
    ///   Required changes in GenerateFromLowess:
    ///     1. Add parallel decoyRtMin/Max/HasWindow arrays (same size as targetRtMin arrays,
    ///        indexed by decoy position: decoyIdx = i - firstDecoyIdx).
    ///     2. Pre-computation loop: for each decoy with RetentionTime, compute its own
    ///        predictedRt / hw / rtMin / rtMax via LOWESS, store in decoy arrays.
    ///     3. Fill loop (decoy path): read from decoyRtMin/Max[decoyIdx] if hasWindow,
    ///        otherwise fall back to targetRtMin/Max[pairedIdx].
    ///
    /// Test structure:
    ///   Layer 1 — Decoy uses its own RT window              [RED until fix]
    ///   Layer 2 — Fallback when decoy has no RT             [GREEN already — preserved by fix]
    ///   Layer 3 — Target window is unchanged by fix         [GREEN already — regression guard]
    ///   Layer 4 — Paired target and decoy windows differ    [RED until fix]
    ///
    /// Placement: mzLib/Test/Dia/DecoyRtWindowTests.cs
    /// </summary>
    [TestFixture]
    public class DecoyRtWindowTests
    {
        // ─────────────────────────────────────────────────────────────────────
        //  LOWESS model factory
        //  Build a minimal valid LowessRtModel on a perfectly linear dataset:
        //    obs = slope × lib + intercept
        //  With 50 points, Fit() always succeeds (minimum is 50).
        //  Because the data is perfectly linear with zero noise:
        //    ToMinutes(lib) = slope × lib + intercept   (exact, no interpolation error)
        //    GetLocalSigma(lib) ≈ 0                     (no residuals)
        //    hw = Math.Min(~0, localSigmaCap) × 2 ≈ 0  → window is essentially a point
        //
        //  To get a non-degenerate hw we use localSigmaCap directly (the cap dominates
        //  when sigma < cap, which it always is for perfect data).
        //  Default localSigmaCap in GenerateFromLowess is 0.5 min → hw = 1.0 min.
        // ─────────────────────────────────────────────────────────────────────

        private const double Slope = 0.5;
        private const double Intercept = 5.0;
        private const double LocalSigmaCap = 0.5;  // default in GenerateFromLowess

        /// <summary>
        /// Fits a perfect linear LOWESS on 50 evenly-spaced anchors:
        ///   lib ∈ [10, 109],  obs = 0.5 × lib + 5.0
        /// ToMinutes(lib) is exact to within floating-point for this data.
        /// </summary>
        private static LowessRtModel BuildLinearLowess()
        {
            int n = 50;
            var libRts = new double[n];
            var obsRts = new double[n];
            for (int i = 0; i < n; i++)
            {
                libRts[i] = 10.0 + i * 2.0;          // 10, 12, 14, … 108
                obsRts[i] = Slope * libRts[i] + Intercept;
            }
            var model = LowessRtModel.Fit(libRts, obsRts, bandwidth: 0.3, enforceMonotonic: true);
            Assert.That(model, Is.Not.Null, "LowessRtModel.Fit must succeed with 50 linear anchors");
            return model;
        }

        /// <summary>
        /// Expected window center for a given library RT through the linear LOWESS model.
        /// </summary>
        private static double ExpectedCenter(double libRt)
            => Slope * libRt + Intercept;

        /// <summary>
        /// Expected half-width: sigma ≈ 0 (perfect data), so cap dominates.
        /// hw = Math.Min(sigma, localSigmaCap) × 2 ≈ localSigmaCap × 2
        /// We use a tolerance of 0.2 min in assertions to absorb tiny LOWESS residuals.
        /// </summary>
        private const double HwApprox = LocalSigmaCap * 2.0;   // ≈ 1.0 min
        private const double WindowTolerance = 0.2;            // min; absorbs tiny residuals

        // ─────────────────────────────────────────────────────────────────────
        //  Scan index factory — one window covering [540, 560] Da, RT [0, 120] min
        // ─────────────────────────────────────────────────────────────────────

        private static DiaScanIndex BuildScanIndex()
        {
            double windowCenter = 550.0;
            double windowWidth = 20.0;

            var scans = new[]
            {
                new MsDataScan(
                    massSpectrum:    new MzSpectrum(new double[] { 100.0 }, new double[] { 1000.0 }, false),
                    oneBasedScanNumber: 1,
                    msnOrder:        2,
                    isCentroid:      true,
                    polarity:        Polarity.Positive,
                    retentionTime:   0.0,
                    scanWindowRange: new MzRange(100, 2000),
                    scanFilter:      "FTMS",
                    mzAnalyzer:      MZAnalyzerType.Orbitrap,
                    totalIonCurrent: 1000.0,
                    injectionTime:   20.0,
                    noiseData:       null,
                    nativeId:        "scan=1",
                    isolationMZ:     windowCenter,
                    isolationWidth:  windowWidth,
                    dissociationType: DissociationType.HCD),

                new MsDataScan(
                    massSpectrum:    new MzSpectrum(new double[] { 100.0 }, new double[] { 1000.0 }, false),
                    oneBasedScanNumber: 2,
                    msnOrder:        2,
                    isCentroid:      true,
                    polarity:        Polarity.Positive,
                    retentionTime:   120.0,
                    scanWindowRange: new MzRange(100, 2000),
                    scanFilter:      "FTMS",
                    mzAnalyzer:      MZAnalyzerType.Orbitrap,
                    totalIonCurrent: 1000.0,
                    injectionTime:   20.0,
                    noiseData:       null,
                    nativeId:        "scan=2",
                    isolationMZ:     windowCenter,
                    isolationWidth:  windowWidth,
                    dissociationType: DissociationType.HCD),
            };

            return DiaScanIndexBuilder.Build(scans);
        }

        private static LibraryPrecursorInput MakePrecursor(
            double precursorMz, double? retentionTime, bool isDecoy, string sequence)
        {
            return new LibraryPrecursorInput(
                sequence: sequence,
                precursorMz: precursorMz,
                chargeState: 2,
                retentionTime: retentionTime,
                isDecoy: isDecoy,
                fragmentMzs: new float[] { 100f, 200f, 300f },
                fragmentIntensities: new float[] { 900f, 500f, 100f },
                irtValue: null);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Helpers to pull a named group from GenerationResult
        // ─────────────────────────────────────────────────────────────────────

        private static DiaLibraryQueryGenerator.PrecursorQueryGroup GetGroup(
            DiaLibraryQueryGenerator.GenerationResult result, int inputIndex)
        {
            foreach (var g in result.PrecursorGroups)
                if (g.InputIndex == inputIndex) return g;
            throw new InvalidOperationException($"No group with InputIndex={inputIndex}");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 1: Decoy uses its own RT window
        //  RED until fix — currently decoy gets target's window
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("RED until fix: decoy RT window must be centered on the decoy's own RetentionTime, not the target's")]
        public void GenerateFromLowess_DecoyWithOwnRt_WindowCenteredOnDecoyRt()
        {
            // Target libRT=30  → predicted obs RT = 0.5×30 + 5 = 20.0 min
            // Decoy  libRT=70  → predicted obs RT = 0.5×70 + 5 = 40.0 min
            // After fix: decoy window is centered at ~40, not ~20.
            // Before fix: decoy window is centered at ~20 (same as target).

            double targetLibRt = 30.0;
            double decoyLibRt = 70.0;
            double expectedDecoyCenter = ExpectedCenter(decoyLibRt);  // 40.0
            double expectedTargetCenter = ExpectedCenter(targetLibRt); // 20.0

            Assert.That(expectedDecoyCenter, Is.Not.EqualTo(expectedTargetCenter).Within(1.0),
                "Test setup: target and decoy must have meaningfully different predicted RTs");

            using var scanIndex = BuildScanIndex();
            var lowess = BuildLinearLowess();

            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.0, targetLibRt, isDecoy: false, sequence: "TARGETK"),
                MakePrecursor(550.0, decoyLibRt,  isDecoy: true,  sequence: "TEGRATEK"),
            };

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, lowess,
                localSigmaCap: LocalSigmaCap);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(2),
                "Both target and decoy must produce a group");

            var targetGroup = GetGroup(result, inputIndex: 0);
            var decoyGroup = GetGroup(result, inputIndex: 1);

            double decoyCenter = (decoyGroup.RtMin + decoyGroup.RtMax) / 2.0;
            double targetCenter = (targetGroup.RtMin + targetGroup.RtMax) / 2.0;

            // Core assertion: decoy window is centered on its own RT, not the target's
            Assert.That(decoyCenter, Is.EqualTo(expectedDecoyCenter).Within(WindowTolerance),
                $"Decoy window center ({decoyCenter:F3}) must be near predicted decoy RT " +
                $"({expectedDecoyCenter:F3}), not target RT ({expectedTargetCenter:F3}). " +
                "Fails today because GenerateFromLowess uses the paired target's window.");

            // Also verify the target window is unchanged
            Assert.That(targetCenter, Is.EqualTo(expectedTargetCenter).Within(WindowTolerance),
                $"Target window center must remain at {expectedTargetCenter:F3}");
        }

        [Test]
        [Description("RED until fix: decoy window center must differ from target window center when RTs differ")]
        public void GenerateFromLowess_DecoyWithDifferentRt_WindowsDiffer()
        {
            // Simpler version of the above: just assert the two windows are different.
            // This is the minimal regression guard — if this fails, the feature is broken.
            double targetLibRt = 30.0;
            double decoyLibRt = 70.0;

            using var scanIndex = BuildScanIndex();
            var lowess = BuildLinearLowess();

            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.0, targetLibRt, isDecoy: false, sequence: "TARGETK"),
                MakePrecursor(550.0, decoyLibRt,  isDecoy: true,  sequence: "TEGRATEK"),
            };

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, lowess,
                localSigmaCap: LocalSigmaCap);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(2));

            var targetGroup = GetGroup(result, inputIndex: 0);
            var decoyGroup = GetGroup(result, inputIndex: 1);

            // Decoy and target windows must not be identical when their RTs differ by 20+ min
            Assert.That(decoyGroup.RtMin, Is.Not.EqualTo(targetGroup.RtMin).Within(WindowTolerance),
                "Decoy RtMin must differ from target RtMin after fix. " +
                "Currently identical because decoy copies the paired target's window.");
            Assert.That(decoyGroup.RtMax, Is.Not.EqualTo(targetGroup.RtMax).Within(WindowTolerance),
                "Decoy RtMax must differ from target RtMax after fix.");
        }

        [Test]
        [Description("RED until fix: multiple decoys each get their own RT window, not all the same target window")]
        public void GenerateFromLowess_MultipleDecoys_EachGetsOwnWindow()
        {
            // Three targets at libRT 20, 50, 80 → predicted obs RT 15, 30, 45
            // Three decoys paired to them at libRT 80, 50, 20 (reversed)
            // After fix: each decoy window is centered on its own predicted RT.
            // Before fix: decoy[i] gets the window of target[i], so:
            //   decoy[0] (libRT=80) gets target[0] window (center ~15) — wrong
            //   decoy[1] (libRT=50) gets target[1] window (center ~30) — coincidentally correct
            //   decoy[2] (libRT=20) gets target[2] window (center ~45) — wrong

            double[] targetLibRts = { 20.0, 50.0, 80.0 };
            double[] decoyLibRts = { 80.0, 50.0, 20.0 }; // reversed

            using var scanIndex = BuildScanIndex();
            var lowess = BuildLinearLowess();

            var precursors = new List<LibraryPrecursorInput>();
            for (int i = 0; i < 3; i++)
                precursors.Add(MakePrecursor(550.0, targetLibRts[i], false, $"TARGET{i}K"));
            for (int i = 0; i < 3; i++)
                precursors.Add(MakePrecursor(550.0, decoyLibRts[i], true, $"YTEGRATEK{i}"));

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, lowess,
                localSigmaCap: LocalSigmaCap);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(6), "All 6 precursors must produce groups");

            // Check each decoy window is centered near its own expected RT
            for (int i = 0; i < 3; i++)
            {
                int decoyInputIdx = 3 + i; // decoys start at index 3
                var decoyGroup = GetGroup(result, decoyInputIdx);
                double decoyCenter = (decoyGroup.RtMin + decoyGroup.RtMax) / 2.0;
                double expectedCenter = ExpectedCenter(decoyLibRts[i]);

                Assert.That(decoyCenter, Is.EqualTo(expectedCenter).Within(WindowTolerance),
                    $"Decoy[{i}] (libRT={decoyLibRts[i]}) window center ({decoyCenter:F3}) " +
                    $"must be near its own predicted RT ({expectedCenter:F3}). " +
                    $"Fails today because it copies the paired target[{i}] window " +
                    $"(center ~{ExpectedCenter(targetLibRts[i]):F3}).");
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 2: Fallback when decoy has no RetentionTime
        //  GREEN already — this is preserved by the fix
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("GREEN: Decoy with no RetentionTime must fall back to paired target's window")]
        public void GenerateFromLowess_DecoyWithNoRt_FallsBackToTargetWindow()
        {
            // Target has libRT=40 → predicted obs RT = 0.5×40 + 5 = 25.0
            // Decoy has no RT → must fall back to target window (center ~25)
            double targetLibRt = 40.0;
            double expectedTargetCenter = ExpectedCenter(targetLibRt); // 25.0

            using var scanIndex = BuildScanIndex();
            var lowess = BuildLinearLowess();

            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.0, targetLibRt, isDecoy: false, sequence: "TARGETK"),
                MakePrecursor(550.0, null,        isDecoy: true,  sequence: "TEGRATEK"),
            };

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, lowess,
                localSigmaCap: LocalSigmaCap);

            // Decoy with no RT: falls back to target window
            // (currently skipped with skippedNoFragments++ in the no-RT path;
            //  if the fix correctly falls back, the decoy group must appear)
            var targetGroup = GetGroup(result, inputIndex: 0);
            double targetCenter = (targetGroup.RtMin + targetGroup.RtMax) / 2.0;

            Assert.That(targetCenter, Is.EqualTo(expectedTargetCenter).Within(WindowTolerance),
                "Target window must be correct regardless of fix");

            // Decoy group: must be present (fallback to target window, not skipped)
            Assert.That(result.PrecursorGroups.Any(g => g.InputIndex == 1), Is.True,
                "Decoy with no RT must fall back to paired target window and still appear in results");

            var decoyGroup = GetGroup(result, inputIndex: 1);
            double decoyCenter = (decoyGroup.RtMin + decoyGroup.RtMax) / 2.0;

            Assert.That(decoyCenter, Is.EqualTo(expectedTargetCenter).Within(WindowTolerance),
                "Decoy with no RT must fall back to the target window center");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 3: Target windows are unchanged by the fix
        //  GREEN already — regression guard
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("GREEN: Target RT windows must be unaffected by the decoy window fix")]
        public void GenerateFromLowess_TargetWindows_UnchangedByFix()
        {
            // Three targets with well-separated RTs. Their windows must be exactly as
            // computed from LOWESS regardless of what decoys do.
            double[] targetLibRts = { 20.0, 50.0, 80.0 };

            using var scanIndex = BuildScanIndex();
            var lowess = BuildLinearLowess();

            var precursors = new List<LibraryPrecursorInput>();
            for (int i = 0; i < 3; i++)
                precursors.Add(MakePrecursor(550.0, targetLibRts[i], false, $"TARGET{i}K"));
            // Decoys with own RTs (reversed) to force a behavioral difference
            double[] decoyLibRts = { 80.0, 50.0, 20.0 };
            for (int i = 0; i < 3; i++)
                precursors.Add(MakePrecursor(550.0, decoyLibRts[i], true, $"YTEGRATEK{i}"));

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, lowess,
                localSigmaCap: LocalSigmaCap);

            for (int i = 0; i < 3; i++)
            {
                var targetGroup = GetGroup(result, inputIndex: i);
                double targetCenter = (targetGroup.RtMin + targetGroup.RtMax) / 2.0;
                double expectedCenter = ExpectedCenter(targetLibRts[i]);

                Assert.That(targetCenter, Is.EqualTo(expectedCenter).Within(WindowTolerance),
                    $"Target[{i}] window center must be {expectedCenter:F3} regardless of decoy fix");
            }
        }

        [Test]
        [Description("GREEN: Decoy with same RT as its target gets same window (degenerate case)")]
        public void GenerateFromLowess_DecoyWithSameRtAsTarget_SameWindow()
        {
            // If decoy libRT == target libRT, the window should be identical.
            // This is a degenerate case (decoy and target happen to have the same library RT)
            // but the fix must handle it without error.
            double sharedLibRt = 40.0;
            double expectedCenter = ExpectedCenter(sharedLibRt); // 25.0

            using var scanIndex = BuildScanIndex();
            var lowess = BuildLinearLowess();

            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.0, sharedLibRt, isDecoy: false, sequence: "TARGETK"),
                MakePrecursor(550.0, sharedLibRt, isDecoy: true,  sequence: "TEGRATEK"),
            };

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, lowess,
                localSigmaCap: LocalSigmaCap);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(2));

            var targetGroup = GetGroup(result, inputIndex: 0);
            var decoyGroup = GetGroup(result, inputIndex: 1);

            double targetCenter = (targetGroup.RtMin + targetGroup.RtMax) / 2.0;
            double decoyCenter = (decoyGroup.RtMin + decoyGroup.RtMax) / 2.0;

            Assert.That(targetCenter, Is.EqualTo(expectedCenter).Within(WindowTolerance));
            Assert.That(decoyCenter, Is.EqualTo(expectedCenter).Within(WindowTolerance),
                "Decoy with same libRT as target must produce same window center");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 4: Decoy window width is independently computed (not inherited)
        //  RED until fix — verifies the hw computation runs for the decoy's own RT
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        [Description("RED until fix: decoy window width must be independently computed from decoy's local sigma")]
        public void GenerateFromLowess_DecoyWindow_WidthDerivedFromDecoyLocalSigma()
        {
            // With perfect linear LOWESS data, sigma ≈ 0 at every point, so
            // hw = min(sigma, cap) × 2 ≈ 0 and windows are near-zero width.
            // To get a non-degenerate window we need sigma > 0.
            //
            // Strategy: use a LOWESS model with deliberate noise so sigma is large
            // enough that cap does NOT dominate (sigma > cap → hw = cap × 2 = 1.0 min).
            // We add ±2 min Gaussian-like noise to the obs RTs, which gives sigma > 0.5.
            //
            // What the test actually pins: both target and decoy independently run the
            // same sigma/hw computation and produce non-zero-width windows.
            // Before fix: decoy inherits the target's rtMin/rtMax values, which may
            // have been computed from a different (possibly zero-sigma) local region.

            var rng = new Random(42);
            int n = 100; // well above 50-anchor minimum
            var libRts = new double[n];
            var obsRts = new double[n];
            for (int i = 0; i < n; i++)
            {
                libRts[i] = 10.0 + i * 1.0;
                // Add ±2 min noise so residuals are large and sigma > localSigmaCap
                double noise = (rng.NextDouble() * 4.0) - 2.0;
                obsRts[i] = Slope * libRts[i] + Intercept + noise;
            }
            var noisyLowess = LowessRtModel.Fit(libRts, obsRts, bandwidth: 0.3, enforceMonotonic: false);
            Assert.That(noisyLowess, Is.Not.Null);

            // Confirm sigma exceeds our cap so windows will be non-degenerate
            double sigmaAtTarget = noisyLowess.GetLocalSigma(30.0);
            double sigmaAtDecoy = noisyLowess.GetLocalSigma(70.0);
            Assert.That(sigmaAtTarget, Is.GreaterThan(LocalSigmaCap),
                "Test setup: target sigma must exceed cap so hw > 0");
            Assert.That(sigmaAtDecoy, Is.GreaterThan(LocalSigmaCap),
                "Test setup: decoy sigma must exceed cap so hw > 0");

            double targetLibRt = 30.0;
            double decoyLibRt = 70.0;

            using var scanIndex = BuildScanIndex();

            var precursors = new List<LibraryPrecursorInput>
            {
                MakePrecursor(550.0, targetLibRt, isDecoy: false, sequence: "TARGETK"),
                MakePrecursor(550.0, decoyLibRt,  isDecoy: true,  sequence: "TEGRATEK"),
            };

            var result = DiaLibraryQueryGenerator.GenerateFromLowess(
                precursors, scanIndex, ppmTolerance: 10f, noisyLowess,
                localSigmaCap: LocalSigmaCap);

            Assert.That(result.PrecursorGroups.Length, Is.EqualTo(2));

            var targetGroup = GetGroup(result, inputIndex: 0);
            var decoyGroup = GetGroup(result, inputIndex: 1);

            double targetWidth = targetGroup.RtMax - targetGroup.RtMin;
            double decoyWidth = decoyGroup.RtMax - decoyGroup.RtMin;

            // Both windows must be non-degenerate (sigma > cap → hw = cap × 2 = 1.0 min)
            double expectedWidth = LocalSigmaCap * 2.0 * 2.0; // hw = cap×2, width = 2×hw
            Assert.That(targetWidth, Is.EqualTo(expectedWidth).Within(WindowTolerance),
                $"Target window width should be ~{expectedWidth:F2} min (cap-dominated)");
            Assert.That(decoyWidth, Is.EqualTo(expectedWidth).Within(WindowTolerance),
                $"Decoy window width should be ~{expectedWidth:F2} min independently computed. " +
                "Fails today if decoy inherits the target's pre-computed window bounds.");

            // And the centers must differ (decoy at its own RT, not target's)
            double targetCenter = (targetGroup.RtMin + targetGroup.RtMax) / 2.0;
            double decoyCenter = (decoyGroup.RtMin + decoyGroup.RtMax) / 2.0;
            Assert.That(decoyCenter, Is.Not.EqualTo(targetCenter).Within(WindowTolerance),
                "Decoy window center must differ from target window center");
        }
    }
}