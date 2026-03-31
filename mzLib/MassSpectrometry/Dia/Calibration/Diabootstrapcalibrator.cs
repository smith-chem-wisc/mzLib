// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: mzLib/MassSpectrometry/Dia/Calibration/DiaBootstrapCalibrator.cs
//
// KEY DESIGN DECISIONS (v5)
//
// 1. PEAK WIDTH FLOOR
//    The half-width floor is derived from the data itself: after each step, the
//    median FWHM of confident target peaks is measured from the TIC array.
//    HalfWidth = max(residualSigma × 2,  medianFWHM × 1.5)
//    This prevents the window from tightening below a physically meaningful
//    minimum regardless of how small σ gets.
//
// 2. SLIDING ANCHOR WINDOW
//    The RT line is fitted from anchors pooled across the last 3 confident steps
//    (those with R² ≥ MinR2ForSlopeTrust).  This gives the line a long RT
//    baseline to fit against, making the slope estimate much more stable.
//
// 3. SLOPE RESCUE
//    If the current step's R² < MinR2ForSlopeTrust, the slope from the pooled
//    prior anchors is kept; only the intercept is re-estimated from the current
//    step's anchors.  This keeps the line tracking the curve even in sparse
//    early-gradient regions where the library RT axis is poorly correlated.

using System;
using System.Buffers;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia.Calibration
{
    public static class DiaBootstrapCalibrator
    {
        // ── Tuning ──────────────────────────────────────────────────────────
        public const float SliceFraction = 0.20f;
        public const float StepPct = 0.05f;
        public const double DefaultLowessBandwidth = 0.3;
        public const int MaxIterations = 15;
        public const float ConvergenceEpsilon = 0.002f;
        public const int MinTargetsRequired = 30;
        public const int MinDecoysRequired = 10;
        public const float SnrGateThreshold = 2.0f;
        public const int MinGatedFragments = 3;
        public const int FeatureCount = 15;

        /// <summary>
        /// R² below this is considered a poor fit — slope is rescued from prior anchors.
        /// </summary>
        public const float MinR2ForSlopeTrust = 0.20f;

        /// <summary>
        /// FWHM multiplier for the half-width floor.
        /// HalfWidth ≥ medianFWHM × FwhmFloorMultiplier.
        /// </summary>
        public const float FwhmFloorMultiplier = 1.5f;

        /// <summary>Number of prior confident steps whose anchors are pooled for the RT line.</summary>
        public const int SlidingAnchorWindow = 3;

        public const int FeatApexDotProduct = 0;
        public const int FeatMeanFragCorr = 1;
        public const int FeatMinFragCorr = 2;
        public const int FeatLog2SNR = 3;
        public const int FeatApexToMeanRatio = 4;
        public const int FeatBoundarySignalRatio = 5;
        public const int FeatMeanSignalRatioDev = 6;
        public const int FeatMeanMassErrorPpm = 7;
        public const int FeatFragDetectionRate = 8;
        public const int FeatBestFragCorrSum = 9;
        public const int FeatMedianFragRefCorr = 10;
        public const int FeatMs1Ms2Correlation = 11;
        public const int FeatCandidateCount = 12;
        public const int FeatRtDeviationMinutes = 13;
        public const int FeatRtDeviationSquared = 14;

        public static readonly string[] FeatureNames =
        {
            "ApexDotProduct", "MeanFragCorr", "MinFragCorr",
            "Log2SNR", "ApexToMeanRatio", "BoundarySignalRatio",
            "MeanSignalRatioDev", "MeanMassErrorPpm", "FragDetRate",
            "BestFragCorrSum", "MedianFragRefCorr", "Ms1Ms2Corr",
            "CandidateCount", "RtDeviationMin", "RtDeviationSq"
        };

        // ══════════════════════════════════════════════════════════════════
        // Bidirectional bootstrap — production entry point
        // ══════════════════════════════════════════════════════════════════

        /// <summary>
        /// Runs a full bidirectional sliding-window bootstrap calibration.
        ///
        /// Algorithm:
        ///   1. Seed at the middle <see cref="SliceFraction"/> of the run (no RT restriction).
        ///   2. Left arm: shift left by <paramref name="stepPct"/> per step until the left
        ///      edge, each step using warm-start weights and a sliding anchor deque.
        ///   3. Right arm: same, shifting right.
        ///   4. Pool all anchors from seed + both arms.
        ///   5. Fit a global OLS line (baseline / diagnostics).
        ///   6. Fit a LOWESS model (primary calibration output).
        ///
        /// Returns a <see cref="BootstrapCalibrationResult"/> whose
        /// <see cref="BootstrapCalibrationResult.LowessModel"/> is the primary artifact.
        /// Use it to generate per-precursor RT windows for the full-run search:
        /// <code>
        ///   double predictedRt = result.LowessModel.ToMinutes(precursor.LibraryRt);
        ///   double hw = Math.Min(result.LowessModel.GetLocalSigma(precursor.LibraryRt), 1.5) * 2.0;
        /// </code>
        /// </summary>
        public static BootstrapCalibrationResult RunBidirectional(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            float ppmTolerance = 20f,
            float stepPct = StepPct,
            double lowessBandwidth = DefaultLowessBandwidth,
            Action<string>? progressReporter = null,
            float seedFwhmMultiplier = 4.0f,
            bool armsEnabled = true)
        {
            void Log(string msg) => progressReporter?.Invoke(msg);

            // ── Seed pass 1 — cold start, no RT prediction ───────────────────
            (float seedLo, float seedHi) = ComputeSliceBounds(scanIndex);
            Log($"[Bootstrap] Seed pass 1: [{seedLo:F3}, {seedHi:F3}] min");

            BootstrapSliceResult? seed1 = RunSlice(
                precursors, scanIndex, seedLo, seedHi, ppmTolerance,
                warmStartWeights: null,
                progressReporter: progressReporter,
                rtPrediction: null,
                pooledPriorAnchors: null);

            if (seed1 == null)
            {
                Log("[Bootstrap] Seed pass 1 failed — aborting.");
                return new BootstrapCalibrationResult(
                    null, null, new List<BootstrapAnchor>(),
                    null, new List<BootstrapSliceResult>(), new List<BootstrapSliceResult>());
            }

            Log($"[Bootstrap] Seed pass 1 anchors: {seed1.Anchors.Count}");

            // ── Seed pass 2 — same window, now with RT prediction ────────────
            // Pass 1 gives us a rough RT line. Pass 2 re-extracts the same slice
            // with rtPrediction populated so:
            //   (a) decoys get ±1 min jitter → realistic RT deviation signal
            //   (b) RtDeviationMinutes and RtDeviationSquared are meaningful features
            // The window itself is unchanged — same seedLo/seedHi, same hw from the
            // RT line. No tightening, no floor games.
            BootstrapSliceResult seed;
            if (seed1.RtLine != null && seed1.Anchors.Count >= 5)
            {
                BootstrapSliceResult? seed2 = RunSlice(
                    precursors, scanIndex, seedLo, seedHi, ppmTolerance,
                    warmStartWeights: seed1.Model.Weights,
                    progressReporter: progressReporter,
                    rtPrediction: seed1.RtLine,
                    pooledPriorAnchors: null);

                if (seed2 != null && seed2.Anchors.Count >= seed1.Anchors.Count)
                {
                    Log($"[Bootstrap] Seed pass 2 anchors: {seed2.Anchors.Count}  — using pass 2");
                    seed = seed2;
                }
                else
                {
                    Log($"[Bootstrap] Seed pass 2 anchors: {seed2?.Anchors.Count ?? 0}  — falling back to pass 1");
                    seed = seed1;
                }
            }
            else
            {
                Log($"[Bootstrap] Seed pass 2 skipped (insufficient anchors or no RT line)");
                seed = seed1;
            }

            // ── Arms ────────────────────────────────────────────────────────
            List<BootstrapSliceResult> leftSteps, rightSteps;
            if (armsEnabled)
            {
                Log("[Bootstrap] Starting left arm...");
                leftSteps = RunArm(seed, precursors, scanIndex, ppmTolerance,
                                   stepPct, goRight: false, progressReporter);

                Log("[Bootstrap] Starting right arm...");
                rightSteps = RunArm(seed, precursors, scanIndex, ppmTolerance,
                                    stepPct, goRight: true, progressReporter);
            }
            else
            {
                Log("[Bootstrap] Arms disabled — seed only.");
                leftSteps = new List<BootstrapSliceResult>();
                rightSteps = new List<BootstrapSliceResult>();
            }

            // ── Pool all anchors ─────────────────────────────────────────────
            int cap = seed.Anchors.Count
                    + leftSteps.Sum(s => s.Anchors.Count)
                    + rightSteps.Sum(s => s.Anchors.Count);

            var allAnchors = new List<BootstrapAnchor>(cap);
            allAnchors.AddRange(seed.Anchors);
            foreach (var s in leftSteps) allAnchors.AddRange(s.Anchors);
            foreach (var s in rightSteps) allAnchors.AddRange(s.Anchors);

            Log($"[Bootstrap] Total anchors: {allAnchors.Count}" +
                $"  (seed={seed.Anchors.Count}" +
                $"  left={leftSteps.Sum(s => s.Anchors.Count)}" +
                $"  right={rightSteps.Sum(s => s.Anchors.Count)})");

            // ── Global OLS line ──────────────────────────────────────────────
            BootstrapRtLine? globalOls = BootstrapRtLine.Fit(allAnchors, fwhmFloor: 0f);
            if (globalOls != null)
                Log($"[Bootstrap] OLS global: {globalOls}");
            else
                Log("[Bootstrap] OLS fit failed (insufficient anchors).");

            // ── LOWESS model ─────────────────────────────────────────────────
            LowessRtModel? lowess = FitLowess(allAnchors, lowessBandwidth, progressReporter);

            return new BootstrapCalibrationResult(
                lowess, globalOls, allAnchors, seed, leftSteps, rightSteps);
        }
        // ══════════════════════════════════════════════════════════════════
        // Arm runner
        // ══════════════════════════════════════════════════════════════════

        private static List<BootstrapSliceResult> RunArm(
            BootstrapSliceResult seed,
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            float ppmTolerance,
            float stepPct,
            bool goRight,
            Action<string>? progressReporter)
        {
            string armName = goRight ? "Right" : "Left";
            void Log(string msg) => progressReporter?.Invoke(msg);

            var results = new List<BootstrapSliceResult>();
            var deque = new Queue<BootstrapSliceResult>();

            if (seed.Anchors.Count >= 5)
                deque.Enqueue(seed);

            BootstrapSliceResult prev = seed;
            int step = 1;

            while (true)
            {
                var next = goRight
                    ? prev.ShiftRight(scanIndex, stepPct)
                    : prev.ShiftLeft(scanIndex, stepPct);

                if (next == null) break;

                float lo = next.Value.Lo, hi = next.Value.Hi;
                Log($"[Bootstrap] {armName} step {step}: [{lo:F3}, {hi:F3}] min");

                List<BootstrapAnchor>? pooled = null;
                if (deque.Count > 0)
                {
                    pooled = new List<BootstrapAnchor>(deque.Sum(s => s.Anchors.Count));
                    foreach (var s in deque) pooled.AddRange(s.Anchors);
                }

                // ── Pass 1: extract with inherited RT line ───────────────────
                var pass1 = RunSlice(
                    precursors, scanIndex, lo, hi, ppmTolerance,
                    warmStartWeights: prev.Model.Weights,
                    progressReporter: progressReporter,
                    rtPrediction: prev.RtLine,
                    pooledPriorAnchors: pooled);

                if (pass1 == null)
                {
                    Log($"[Bootstrap] {armName} step {step} pass 1 produced no model — stopping arm.");
                    break;
                }

                // ── Pass 2: refit local RT line from pass 1 anchors, re-extract
                // This gives decoys a locally-corrected jitter target and gives the
                // LDA a clean RT deviation signal against the local gradient shape.
                BootstrapSliceResult result = pass1;
                if (pass1.RtLine != null && pass1.Anchors.Count >= 5)
                {
                    var pass2 = RunSlice(
                        precursors, scanIndex, lo, hi, ppmTolerance,
                        warmStartWeights: pass1.Model.Weights,
                        progressReporter: progressReporter,
                        rtPrediction: pass1.RtLine,
                        pooledPriorAnchors: pooled);

                    if (pass2 != null && pass2.Anchors.Count >= pass1.Anchors.Count)
                        result = pass2;
                }

                results.Add(result);
                prev = result;

                if (result.Anchors.Count >= 5)
                {
                    deque.Enqueue(result);
                    if (deque.Count > SlidingAnchorWindow)
                        deque.Dequeue();
                }

                step++;
            }

            return results;
        }

        // ══════════════════════════════════════════════════════════════════
        // LOWESS fit
        // ══════════════════════════════════════════════════════════════════

        private static LowessRtModel? FitLowess(
            IReadOnlyList<BootstrapAnchor> anchors,
            double bandwidth,
            Action<string>? progressReporter)
        {
            void Log(string msg) => progressReporter?.Invoke(msg);

            var libRts = new List<double>(anchors.Count);
            var obsRts = new List<double>(anchors.Count);
            foreach (var a in anchors)
                if (a.LibraryRt > 0f && a.ObservedApexRt > 0f)
                {
                    libRts.Add(a.LibraryRt);
                    obsRts.Add(a.ObservedApexRt);
                }

            if (libRts.Count < 50)
            {
                Log($"[Bootstrap] LOWESS skipped — only {libRts.Count} valid anchors (need ≥50).");
                return null;
            }

            Log($"[Bootstrap] Fitting LOWESS on {libRts.Count} anchors, bandwidth={bandwidth:F2}...");
            var sw = System.Diagnostics.Stopwatch.StartNew();

            var model = LowessRtModel.Fit(
                libRts.ToArray(), obsRts.ToArray(),
                bandwidth: bandwidth,
                enforceMonotonic: true);

            sw.Stop();
            Log($"[Bootstrap] LOWESS done in {sw.Elapsed.TotalSeconds:F2}s" +
                $"  R²={model?.RSquared:F4}  σ={model?.SigmaMinutes:F4} min");

            return model;
        }

        // ══════════════════════════════════════════════════════════════════
        // Piecewise prefit + LOWESS residual correction
        // ══════════════════════════════════════════════════════════════════

        /// <summary>
        /// Two-stage calibration that corrects the shape of the mapping before
        /// LOWESS is applied.
        ///
        /// Problem with plain LOWESS on bootstrap anchors:
        ///   The steep region of the gradient (20–35 min) has far more anchors
        ///   than the flat ends. LOWESS uses a fixed bandwidth fraction of N,
        ///   so dense regions dominate the local fits and pull the slope flatter
        ///   than the true curve. The result is a systematic tilt: predictions
        ///   are too early on the right side and too late on the left side.
        ///
        /// Solution:
        ///   Stage 1 — Fit a piecewise linear model using evenly-spaced knots
        ///             across the library RT range (not evenly-spaced anchor counts).
        ///             Equal RT intervals means each segment captures the true local
        ///             slope regardless of anchor density.
        ///   Stage 2 — Compute residuals: obsRT - piecewise(libRT). These are small
        ///             (the shape is already correct) and symmetric around zero.
        ///   Stage 3 — Fit LOWESS on (libRT, residuals). LOWESS now only models
        ///             local deviations from an already-correct shape.
        ///   Stage 4 — Resample the combined model (piecewise + residual LOWESS)
        ///             at all anchor library RTs and fit a final LowessRtModel on
        ///             those resampled predictions. This produces a standard
        ///             LowessRtModel that the rest of the pipeline can use unchanged,
        ///             with GetLocalSigma() reflecting the true residual scatter.
        ///
        /// Falls back to plain FitLowess if the piecewise prefit fails.
        /// </summary>
        private static LowessRtModel? FitLowessWithPiecewisePrefit(
            IReadOnlyList<BootstrapAnchor> anchors,
            double bandwidth,
            Action<string>? progressReporter,
            int evenKnots = 25,
            double fullLibRtMin = 0,
            double fullLibRtMax = 170)
        {
            void Log(string msg) => progressReporter?.Invoke(msg);

            // Collect valid anchors
            var libRts = new List<double>(anchors.Count);
            var obsRts = new List<double>(anchors.Count);
            foreach (var a in anchors)
                if (a.LibraryRt > 0f && a.ObservedApexRt > 0f)
                {
                    libRts.Add(a.LibraryRt);
                    obsRts.Add(a.ObservedApexRt);
                }

            if (libRts.Count < 50)
            {
                Log($"[Bootstrap] LOWESS+prefit skipped — only {libRts.Count} valid anchors (need ≥50).");
                return null;
            }

            double[] libArr = libRts.ToArray();
            double[] obsArr = obsRts.ToArray();

            // ── Stage 1: Piecewise linear with evenly-spaced RT knots ────────
            var sw = System.Diagnostics.Stopwatch.StartNew();
            double[]? piecewisePred = FitPiecewiseEvenKnots(
                libArr, obsArr, evenKnots, Log,
                libRtMin: fullLibRtMin, libRtMax: fullLibRtMax);

            if (piecewisePred == null)
            {
                Log("[Bootstrap] Piecewise prefit failed — falling back to plain LOWESS.");
                return FitLowess(anchors, bandwidth, progressReporter);
            }

            sw.Stop();
            Log($"[Bootstrap] Piecewise prefit done in {sw.Elapsed.TotalSeconds:F2}s  " +
                $"knots={evenKnots}");

            // ── Stage 2: Residuals ────────────────────────────────────────────
            var residuals = new double[libArr.Length];
            double ssResPrefit = 0, ssTot = 0, meanObs = 0;
            for (int i = 0; i < obsArr.Length; i++) meanObs += obsArr[i];
            meanObs /= obsArr.Length;
            for (int i = 0; i < libArr.Length; i++)
            {
                residuals[i] = obsArr[i] - piecewisePred[i];
                ssResPrefit += residuals[i] * residuals[i];
                double dev = obsArr[i] - meanObs;
                ssTot += dev * dev;
            }
            double sigmaPrefit = Math.Sqrt(ssResPrefit / libArr.Length);
            double r2Prefit = ssTot > 0 ? 1.0 - ssResPrefit / ssTot : 0;
            Log($"[Bootstrap] Piecewise prefit: R²={r2Prefit:F4}  σ={sigmaPrefit:F4} min  " +
                $"(residuals to LOWESS are much smaller than raw obsRT)");

            // ── Stage 3: LOWESS on residuals ──────────────────────────────────
            Log($"[Bootstrap] Fitting LOWESS on {libArr.Length} residuals, bandwidth={bandwidth:F2}...");
            sw.Restart();

            var residualLowess = LowessRtModel.Fit(
                libArr, residuals,
                bandwidth: bandwidth,
                enforceMonotonic: false);   // residuals need NOT be monotonic

            sw.Stop();
            if (residualLowess == null)
            {
                Log("[Bootstrap] Residual LOWESS failed — falling back to plain LOWESS.");
                return FitLowess(anchors, bandwidth, progressReporter);
            }
            Log($"[Bootstrap] Residual LOWESS done in {sw.Elapsed.TotalSeconds:F2}s  " +
                $"R²={residualLowess.RSquared:F4}  σ={residualLowess.SigmaMinutes:F4} min");

            // ── Stage 4: Resample combined model → final LowessRtModel ───────
            // Evaluate combined = piecewise + residualLowess at every anchor libRT.
            // Fit a final monotonic LowessRtModel on those resampled values so the
            // rest of the pipeline gets a standard LowessRtModel with correct
            // GetLocalSigma() values.
            var combinedObs = new double[libArr.Length];
            for (int i = 0; i < libArr.Length; i++)
                combinedObs[i] = piecewisePred[i] + residualLowess.ToMinutes(libArr[i]);

            Log($"[Bootstrap] Fitting final LOWESS on {libArr.Length} combined predictions...");
            sw.Restart();

            var finalLowess = LowessRtModel.Fit(
                libArr, combinedObs,
                bandwidth: bandwidth,
                enforceMonotonic: true);

            sw.Stop();
            if (finalLowess == null)
            {
                Log("[Bootstrap] Final LOWESS failed — falling back to plain LOWESS.");
                return FitLowess(anchors, bandwidth, progressReporter);
            }

            // The sigma of the final model should be measured against the original
            // observed RTs, not the resampled combined predictions.
            // Report both for diagnostics.
            double ssResFinal = 0;
            for (int i = 0; i < libArr.Length; i++)
            {
                double pred = finalLowess.ToMinutes(libArr[i]);
                double r = obsArr[i] - pred;
                ssResFinal += r * r;
            }
            double sigmaFinal = Math.Sqrt(ssResFinal / libArr.Length);

            Log($"[Bootstrap] Final LOWESS done in {sw.Elapsed.TotalSeconds:F2}s  " +
                $"R²={finalLowess.RSquared:F4}  " +
                $"σ(vs combined)={finalLowess.SigmaMinutes:F4} min  " +
                $"σ(vs original obsRT)={sigmaFinal:F4} min");

            return finalLowess;
        }

        /// <summary>
        /// Fits a piecewise linear mapping using evenly-spaced knots across the
        /// library RT range. This is immune to anchor density imbalance: each
        /// segment spans the same width in library RT space regardless of how
        /// many anchors fall in it.
        ///
        /// Returns predicted observed RTs at each input library RT (same length
        /// as libArr, sorted ascending), or null if fitting fails.
        ///
        /// Each segment is fit by OLS on the anchors within that RT interval.
        /// Segments with fewer than MinAnchorsPerSegment anchors borrow slope
        /// from their neighbors via distance-weighted average.
        /// Monotonicity is enforced at knot boundaries after fitting.
        /// </summary>
        private static double[]? FitPiecewiseEvenKnots(
            double[] libArr,    // sorted ascending
            double[] obsArr,
            int numKnots,
            Action<string> log,
            int minAnchorsPerSegment = 5,
            double libRtMin = 0,
            double libRtMax = 170)
        {
            int n = libArr.Length;
            if (n < 50) return null;

            // Use the FULL library RT range for knot spacing, not just the anchor range.
            // This ensures knots cover the entire gradient even when anchors are absent
            // at the early or late ends (common when those steps yield zero confident IDs).
            double libMin = libRtMin;
            double libMax = libRtMax;
            if (libMax <= libMin) return null;

            int numSegments = Math.Max(2, numKnots - 1);
            double segWidth = (libMax - libMin) / numSegments;

            // Knot positions in library RT space — evenly spaced
            var knotLibRt = new double[numKnots];
            var knotObsRt = new double[numKnots];  // filled after fitting
            var segSlope = new double[numSegments];
            var segIntercept = new double[numSegments];
            var segCount = new int[numSegments];

            for (int k = 0; k < numKnots; k++)
                knotLibRt[k] = libMin + k * segWidth;
            knotLibRt[numKnots - 1] = libMax; // ensure exact right boundary

            // Fit OLS per segment
            for (int s = 0; s < numSegments; s++)
            {
                double lo = knotLibRt[s];
                double hi = knotLibRt[s + 1];

                // Collect anchors in [lo, hi)  (last segment includes hi)
                double sx = 0, sy = 0, sxx = 0, sxy = 0;
                int cnt = 0;
                for (int i = 0; i < n; i++)
                {
                    if (libArr[i] < lo) continue;
                    if (s < numSegments - 1 && libArr[i] >= hi) break;
                    sx += libArr[i]; sy += obsArr[i];
                    sxx += libArr[i] * libArr[i]; sxy += libArr[i] * obsArr[i];
                    cnt++;
                }
                segCount[s] = cnt;

                if (cnt >= 2)
                {
                    double d = cnt * sxx - sx * sx;
                    if (Math.Abs(d) > 1e-12)
                    {
                        segSlope[s] = (cnt * sxy - sx * sy) / d;
                        segIntercept[s] = (sy - segSlope[s] * sx) / cnt;
                        continue;
                    }
                }

                // Too few anchors — use segment midpoint mean with slope=0
                double midLib = (lo + hi) * 0.5;
                double midObs = cnt > 0 ? sy / cnt : (lo + hi) * 0.5; // fallback
                segSlope[s] = 0;
                segIntercept[s] = midObs;
            }

            // Smooth slopes for sparse segments by borrowing from neighbors
            for (int pass = 0; pass < 2; pass++)
            {
                for (int s = 0; s < numSegments; s++)
                {
                    if (segCount[s] >= minAnchorsPerSegment) continue;

                    // Distance-weighted average of neighbor slopes
                    double wSum = 0, slopeSum = 0, interceptSum = 0;
                    for (int nb = 0; nb < numSegments; nb++)
                    {
                        if (nb == s || segCount[nb] < minAnchorsPerSegment) continue;
                        double dist = Math.Abs(nb - s) + 1.0;
                        double w = 1.0 / (dist * dist);
                        wSum += w;
                        slopeSum += w * segSlope[nb];
                        interceptSum += w * segIntercept[nb];
                    }
                    if (wSum > 0)
                    {
                        segSlope[s] = slopeSum / wSum;
                        // Re-anchor intercept to pass through segment midpoint mean
                        double midLib = (knotLibRt[s] + knotLibRt[s + 1]) * 0.5;
                        double midObs = segSlope[s] * midLib + interceptSum / wSum;
                        segIntercept[s] = midObs - segSlope[s] * midLib;
                    }
                }
            }

            // Evaluate at knot boundaries
            for (int k = 0; k < numKnots; k++)
            {
                int s = Math.Min(k, numSegments - 1);
                knotObsRt[k] = segSlope[s] * knotLibRt[k] + segIntercept[s];
            }

            // Enforce monotonicity at knots via pool-adjacent-violators
            for (int k = 1; k < numKnots; k++)
                if (knotObsRt[k] < knotObsRt[k - 1])
                    knotObsRt[k] = knotObsRt[k - 1];

            // Log knot table
            log($"[Bootstrap] Piecewise knots (libRT → predRT  anchors):");
            for (int k = 0; k < numKnots; k++)
            {
                int s = Math.Min(k, numSegments - 1);
                log($"[Bootstrap]   [{k,2}] libRT={knotLibRt[k],7:F2}  predRT={knotObsRt[k]:F3}  " +
                    (k < numSegments ? $"n={segCount[k]}" : ""));
            }

            // Interpolate predictions at every anchor libRT
            var predictions = new double[n];
            for (int i = 0; i < n; i++)
            {
                double x = libArr[i];

                // Find bracketing knots via binary search
                int lo = 0, hi2 = numKnots - 2;
                while (lo < hi2)
                {
                    int mid = (lo + hi2) / 2;
                    if (x < knotLibRt[mid + 1]) hi2 = mid;
                    else lo = mid + 1;
                }

                double leftLib = knotLibRt[lo];
                double rightLib = knotLibRt[lo + 1];
                double leftObs = knotObsRt[lo];
                double rightObs = knotObsRt[lo + 1];

                double frac = (rightLib > leftLib)
                    ? (x - leftLib) / (rightLib - leftLib)
                    : 0.0;
                predictions[i] = leftObs + frac * (rightObs - leftObs);
            }

            return predictions;
        }

        // ══════════════════════════════════════════════════════════════════
        // Convenience entry — middle 20% slice, no warm start
        // ══════════════════════════════════════════════════════════════════

        public static BootstrapSliceResult? Run(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            float ppmTolerance = 20f,
            Action<string>? progressReporter = null)
        {
            (float lo, float hi) = ComputeSliceBounds(scanIndex);
            return RunSlice(precursors, scanIndex, lo, hi, ppmTolerance,
                            null, progressReporter, null, null);
        }

        // ══════════════════════════════════════════════════════════════════
        // Core entry
        // ══════════════════════════════════════════════════════════════════

        /// <param name="pooledPriorAnchors">
        /// Anchors pooled from the last N confident steps.  Used to fit the RT line
        /// when the current step has poor R².  Pass null for the seed step.
        /// </param>
        public static BootstrapSliceResult? RunSlice(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            float sliceLo, float sliceHi,
            float ppmTolerance = 20f,
            float[]? warmStartWeights = null,
            Action<string>? progressReporter = null,
            BootstrapRtLine? rtPrediction = null,
            IReadOnlyList<BootstrapAnchor>? pooledPriorAnchors = null)
        {
            void Log(string msg) => progressReporter?.Invoke($"[Bootstrap] {msg}");

            Log($"RT slice: [{sliceLo:F3}, {sliceHi:F3}] min");
            if (sliceHi <= sliceLo) { Log("Degenerate slice."); return null; }

            // ── Generate queries ────────────────────────────────────────────
            (FragmentQuery[] queries, SliceGroup[] groups) =
                GenerateSliceQueries(precursors, scanIndex, sliceLo, sliceHi,
                                     ppmTolerance, rtPrediction);

            Log($"Generated {queries.Length} queries for {groups.Length} precursors" +
                (rtPrediction != null ? " (RT-filtered)" : ""));
            if (queries.Length == 0) { Log("No queries."); return null; }

            // ── Extract XICs ────────────────────────────────────────────────
            var orchestrator = new DiaExtractionOrchestrator(scanIndex);
            ExtractionResult extr = orchestrator.ExtractAll(queries, -1);
            Log($"Extraction complete: {extr.Results.Length} results");

            // ── Features + apex RTs + per-precursor FWHM ───────────────────
            (float[][] features, bool[] isDecoy,
             float[] libraryRts, float[] observedApexRts,
             float[] peakFwhms,
             int gated, int total) =
                ComputeFeatureVectors(precursors, groups, queries,
                                      extr, scanIndex, ppmTolerance, rtPrediction);

            Log($"SNR gate: {gated}/{total} passed  (threshold={SnrGateThreshold:F1}×)");

            int nT = 0, nD = 0;
            for (int i = 0; i < isDecoy.Length; i++) { if (isDecoy[i]) nD++; else nT++; }
            Log($"Feature vectors: {nT} targets, {nD} decoys");

            if (nT < MinTargetsRequired || nD < MinDecoysRequired)
            {
                Log($"Insufficient data."); return null;
            }

            // ── Fit LDA ─────────────────────────────────────────────────────
            BootstrapLdaModel model = FitIterativeLda(features, isDecoy, warmStartWeights, Log);

            // ── Collect confident anchors for this step ─────────────────────
            var stepAnchors = new List<BootstrapAnchor>();
            var stepFwhms = new List<float>();

            for (int i = 0; i < features.Length; i++)
            {
                float score = model.Score(features[i]);
                if (model.IsTarget(score) && !isDecoy[i] && libraryRts[i] > 0f)
                {
                    stepAnchors.Add(new BootstrapAnchor(
                        libraryRts[i], observedApexRts[i], score));
                    if (!float.IsNaN(peakFwhms[i]) && peakFwhms[i] > 0f)
                        stepFwhms.Add(peakFwhms[i]);
                }
            }

            Log($"Confident anchors: {stepAnchors.Count}");

            // ── Median FWHM → half-width floor ─────────────────────────────
            float medianFwhm = MedianList(stepFwhms);
            float fwhmFloor = medianFwhm > 0f ? medianFwhm * FwhmFloorMultiplier : 0f;

            if (fwhmFloor > 0f)
                Log($"Median peak FWHM: {medianFwhm:F4} min  →  floor: {fwhmFloor:F4} min");

            // ── Fit RT line using pooled anchors ────────────────────────────
            var allAnchors = new List<BootstrapAnchor>(stepAnchors);
            if (pooledPriorAnchors != null) allAnchors.AddRange(pooledPriorAnchors);

            BootstrapRtLine? rtLine = null;
            if (allAnchors.Count >= 5)
            {
                rtLine = BootstrapRtLine.Fit(allAnchors, fwhmFloor);

                // Slope rescue: if pooled line has poor R², borrow slope from prior anchors
                if (rtLine != null && rtLine.RSquared < MinR2ForSlopeTrust
                    && pooledPriorAnchors != null && pooledPriorAnchors.Count >= 5)
                {
                    var priorLine = BootstrapRtLine.Fit(pooledPriorAnchors, fwhmFloor);
                    if (priorLine != null && priorLine.RSquared >= MinR2ForSlopeTrust)
                    {
                        rtLine = BootstrapRtLine.FitWithFixedSlope(
                            stepAnchors.Count >= 5 ? stepAnchors : allAnchors,
                            priorLine.Slope, fwhmFloor);
                        Log($"Slope rescued from prior anchors (R²={priorLine.RSquared:F4})");
                    }
                }

                if (rtLine != null)
                    Log($"RT line: {rtLine}");
                else
                    Log("RT line fit failed.");
            }
            else
            {
                Log("Insufficient anchors for RT line fit (need ≥5).");
            }

            return new BootstrapSliceResult(
                model, rtLine, stepAnchors, sliceLo, sliceHi, medianFwhm);
        }

        // ══════════════════════════════════════════════════════════════════
        // Slice bounds
        // ══════════════════════════════════════════════════════════════════

        public static (float RtLo, float RtHi) ComputeSliceBounds(DiaScanIndex scanIndex)
            => ComputeSliceBoundsAtPercentile(scanIndex,
                   (1f - SliceFraction) / 2f, 1f - (1f - SliceFraction) / 2f);

        public static (float RtLo, float RtHi) ComputeSliceBoundsAtPercentile(
            DiaScanIndex scanIndex, float loPct, float hiPct)
        {
            int n = scanIndex.ScanCount;
            if (n == 0) return (0f, 0f);
            float[] buf = ArrayPool<float>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) buf[i] = scanIndex.GetScanRt(i);
                Array.Sort(buf, 0, n);
                return (buf[(int)(loPct * (n - 1))], buf[(int)(hiPct * (n - 1))]);
            }
            finally { ArrayPool<float>.Shared.Return(buf); }
        }

        // ══════════════════════════════════════════════════════════════════
        // Query generation
        // ══════════════════════════════════════════════════════════════════

        // Maximum RT jitter applied to decoy extraction windows (minutes).
        // Decoys get a deterministic random offset in [-DecoyRtJitterMinutes, +DecoyRtJitterMinutes]
        // so the bootstrap LDA can exploit RT deviation to separate T from D.
        public const float DecoyRtJitterMinutes = 1.0f;

        private static (FragmentQuery[] Queries, SliceGroup[] Groups) GenerateSliceQueries(
            IReadOnlyList<LibraryPrecursorInput> precursors,
            DiaScanIndex scanIndex,
            float sliceLo, float sliceHi, float ppmTol,
            BootstrapRtLine? rtPrediction)
        {
            var qList = new List<FragmentQuery>(precursors.Count * 8);
            var gList = new List<SliceGroup>(precursors.Count);

            for (int pi = 0; pi < precursors.Count; pi++)
            {
                var p = precursors[pi];

                float qLo, qHi;
                if (rtPrediction != null && p.RetentionTime.HasValue)
                {
                    float libRt = (float)p.RetentionTime.Value;
                    float hw = rtPrediction.HalfWidth;
                    float pred = rtPrediction.PredictObservedRt(libRt);
                    qLo = pred - hw;
                    qHi = pred + hw;

                    // Shift decoy windows by a deterministic random RT offset so the
                    // bootstrap LDA can use RT deviation as a discriminating feature.
                    // The offset is derived from the sequence hash — reproducible across runs.
                    if (p.IsDecoy)
                    {
                        int hash = p.Sequence.GetHashCode();
                        // Map hash to [-1, +1] uniformly, then scale to jitter range
                        float jitter = ((hash & 0x7FFFFFFF) / (float)0x7FFFFFFF * 2f - 1f)
                                       * DecoyRtJitterMinutes;
                        qLo += jitter;
                        qHi += jitter;
                    }

                    if (qHi < sliceLo || qLo > sliceHi) continue;
                    qLo = Math.Max(qLo, sliceLo);
                    qHi = Math.Min(qHi, sliceHi);
                }
                else { qLo = sliceLo; qHi = sliceHi; }

                int wid = scanIndex.FindWindowForPrecursorMz(p.PrecursorMz);
                if (wid < 0) continue;
                if (!scanIndex.TryGetScanRangeForWindow(wid, out int wS, out int wC)) continue;

                bool has = false;
                for (int si = wS; si < wS + wC && !has; si++)
                    has = scanIndex.GetScanRt(si) >= qLo && scanIndex.GetScanRt(si) <= qHi;
                if (!has) continue;

                int qo = qList.Count;
                for (int fi = 0; fi < p.FragmentCount; fi++)
                    qList.Add(new FragmentQuery(p.FragmentMzs[fi], ppmTol, qLo, qHi, wid, qList.Count));
                gList.Add(new SliceGroup(pi, qo, p.FragmentCount, wid));
            }
            return (qList.ToArray(), gList.ToArray());
        }

        // ══════════════════════════════════════════════════════════════════
        // Feature vectors + FWHM measurement
        // ══════════════════════════════════════════════════════════════════

        private static (float[][] Features, bool[] IsDecoy,
                        float[] LibraryRts, float[] ObservedApexRts,
                        float[] PeakFwhms,
                        int Gated, int Total)
            ComputeFeatureVectors(
                IReadOnlyList<LibraryPrecursorInput> precursors,
                SliceGroup[] groups,
                FragmentQuery[] queries,
                ExtractionResult extr,
                DiaScanIndex scanIndex,
                float ppmTol,
                BootstrapRtLine? rtPrediction = null)
        {
            var fList = new List<float[]>(groups.Length);
            var dList = new List<bool>(groups.Length);
            var libRts = new List<float>(groups.Length);
            var apexRts = new List<float>(groups.Length);
            var fwhms = new List<float>(groups.Length);
            int gated = 0;

            // Pre-compute per-window precursor count for CandidateCount feature
            var windowCandidateCount = new Dictionary<int, int>(32);
            foreach (var g in groups)
            {
                if (!windowCandidateCount.TryGetValue(g.WindowId, out int cnt))
                    windowCandidateCount[g.WindowId] = 1;
                else
                    windowCandidateCount[g.WindowId] = cnt + 1;
            }

            float[] apexBuf = new float[32];
            float[] expBuf = new float[32];

            foreach (var g in groups)
            {
                var p = precursors[g.PrecursorIndex];
                if (g.QueryCount == 0) continue;

                int refFragIdx = -1, refPts = 0;
                for (int qi = g.QueryOffset; qi < g.QueryOffset + g.QueryCount; qi++)
                {
                    int pts = extr.Results[qi].DataPointCount;
                    if (pts > refPts) { refPts = pts; refFragIdx = qi - g.QueryOffset; }
                }
                if (refFragIdx < 0 || refPts == 0) continue;

                var refResult = extr.Results[g.QueryOffset + refFragIdx];
                int timePoints = refResult.DataPointCount;
                int fc = g.QueryCount;

                float[] mat = ArrayPool<float>.Shared.Rent(timePoints * fc);
                float[] ticBuf = ArrayPool<float>.Shared.Rent(timePoints);
                Array.Clear(mat, 0, timePoints * fc);

                try
                {
                    ReadOnlySpan<float> refRts = new ReadOnlySpan<float>(
                        extr.RtBuffer, refResult.RtBufferOffset, timePoints);

                    for (int qi = g.QueryOffset; qi < g.QueryOffset + g.QueryCount; qi++)
                    {
                        int fi = qi - g.QueryOffset;
                        var res = extr.Results[qi];
                        if (res.DataPointCount == 0) continue;
                        AlignToGrid(refRts,
                            new ReadOnlySpan<float>(extr.RtBuffer, res.RtBufferOffset, res.DataPointCount),
                            new ReadOnlySpan<float>(extr.IntensityBuffer, res.IntensityBufferOffset, res.DataPointCount),
                            mat, fi, fc, timePoints);
                    }

                    Span<float> ticSpan = ticBuf.AsSpan(0, timePoints);
                    float apexTic = 0f; int apexT = 0; float sumTic = 0f;
                    for (int t = 0; t < timePoints; t++)
                    {
                        float tic = 0f; int row = t * fc;
                        for (int f = 0; f < fc; f++) tic += mat[row + f];
                        ticSpan[t] = tic; sumTic += tic;
                        if (tic > apexTic) { apexTic = tic; apexT = t; }
                    }

                    float medianTic = MedianOfSpan(ticSpan, timePoints);
                    float snr = medianTic > 0f ? apexTic / medianTic : 0f;
                    if (snr < SnrGateThreshold) continue;
                    gated++;

                    if (apexBuf.Length < fc) { apexBuf = new float[fc]; expBuf = new float[fc]; }
                    else Array.Clear(apexBuf, 0, fc);

                    int apexRow = apexT * fc;
                    int gatedFc = 0;
                    float noiseFloor = apexTic * 0.05f;
                    for (int f = 0; f < fc; f++)
                    {
                        apexBuf[f] = mat[apexRow + f];
                        expBuf[f] = p.FragmentIntensities[f];
                        if (apexBuf[f] > noiseFloor) gatedFc++;
                    }
                    if (gatedFc < MinGatedFragments) continue;

                    float observedApexRt = refRts[apexT];
                    float fwhm = ComputeFwhm(ticSpan, refRts, apexT, apexTic);

                    float dot = NdpScore(new ReadOnlySpan<float>(apexBuf, 0, fc),
                                         new ReadOnlySpan<float>(expBuf, 0, fc));
                    (float mc, float mnc) = PairwiseXicCorrelation(mat, fc, timePoints);
                    float log2Snr = MathF.Log2(Math.Max(snr, 1f));
                    float meanTic = sumTic / timePoints;
                    float a2m = meanTic > 0f ? apexTic / meanTic : 1f;
                    float bndRatio = apexTic > 0f
                        ? (ticSpan[0] + ticSpan[timePoints - 1]) * 0.5f / apexTic : 1f;
                    float srd = ComputeMeanSignalRatioDev(
                        new ReadOnlySpan<float>(apexBuf, 0, fc),
                        new ReadOnlySpan<float>(expBuf, 0, fc));
                    float massErr = ComputeMassErrorPpm(
                        g, queries, extr, scanIndex, observedApexRt, ppmTol);

                    // ── New features ─────────────────────────────────────────
                    (float bestFragCorrSum, float medianFragRefCorr) =
                        ComputeBestFragFeatures(mat, fc, timePoints);

                    float ms1Ms2Corr = ComputeBootstrapMs1Ms2Correlation(
                        scanIndex, g.WindowId, observedApexRt,
                        extr, g.QueryOffset + refFragIdx, ppmTol);

                    windowCandidateCount.TryGetValue(g.WindowId, out int candCount);

                    var vec = new float[FeatureCount];
                    vec[FeatApexDotProduct] = dot;
                    vec[FeatMeanFragCorr] = float.IsNaN(mc) ? 0f : mc;
                    vec[FeatMinFragCorr] = float.IsNaN(mnc) ? -1f : mnc;
                    vec[FeatLog2SNR] = log2Snr;
                    vec[FeatApexToMeanRatio] = a2m;
                    vec[FeatBoundarySignalRatio] = bndRatio;
                    vec[FeatMeanSignalRatioDev] = float.IsNaN(srd) ? 0f : srd;
                    vec[FeatMeanMassErrorPpm] = massErr;
                    vec[FeatFragDetectionRate] = (float)gatedFc / fc;
                    vec[FeatBestFragCorrSum] = float.IsNaN(bestFragCorrSum) ? 0f : bestFragCorrSum;
                    vec[FeatMedianFragRefCorr] = float.IsNaN(medianFragRefCorr) ? 0f : medianFragRefCorr;
                    vec[FeatMs1Ms2Correlation] = float.IsNaN(ms1Ms2Corr) ? 0f : ms1Ms2Corr;
                    vec[FeatCandidateCount] = candCount;

                    // RT deviation — predicted vs observed apex RT.
                    // Targets cluster near 0; jittered decoys are spread across ±1 min.
                    // Use 0 when no RT prediction is available (cold seed pass).
                    float rtDev = 0f;
                    if (rtPrediction != null && p.RetentionTime.HasValue)
                    {
                        float predicted = rtPrediction.PredictObservedRt((float)p.RetentionTime.Value);
                        rtDev = observedApexRt - predicted;
                    }
                    vec[FeatRtDeviationMinutes] = rtDev;
                    vec[FeatRtDeviationSquared] = rtDev * rtDev;

                    fList.Add(vec);
                    dList.Add(p.IsDecoy);
                    libRts.Add(p.RetentionTime.HasValue ? (float)p.RetentionTime.Value : float.NaN);
                    apexRts.Add(observedApexRt);
                    fwhms.Add(fwhm);
                }
                finally
                {
                    ArrayPool<float>.Shared.Return(mat);
                    ArrayPool<float>.Shared.Return(ticBuf);
                }
            }

            return (fList.ToArray(), dList.ToArray(),
                    libRts.ToArray(), apexRts.ToArray(),
                    fwhms.ToArray(), gated, groups.Length);
        }

        // ══════════════════════════════════════════════════════════════════
        // FWHM computation
        // ══════════════════════════════════════════════════════════════════

        internal static float ComputeFwhm(
            Span<float> tic, ReadOnlySpan<float> rts, int apexT, float apexTic)
        {
            float halfMax = apexTic * 0.5f;
            int n = tic.Length;

            float leftRt = rts[0];
            bool foundLeft = false;
            for (int t = apexT - 1; t >= 0; t--)
            {
                if (tic[t] <= halfMax)
                {
                    float span = tic[t + 1] - tic[t];
                    float frac = span > 0f ? (halfMax - tic[t]) / span : 0.5f;
                    leftRt = rts[t] + frac * (rts[t + 1] - rts[t]);
                    foundLeft = true;
                    break;
                }
            }
            if (!foundLeft) return float.NaN;

            float rightRt = rts[n - 1];
            bool foundRight = false;
            for (int t = apexT + 1; t < n; t++)
            {
                if (tic[t] <= halfMax)
                {
                    float span = tic[t - 1] - tic[t];
                    float frac = span > 0f ? (halfMax - tic[t]) / span : 0.5f;
                    rightRt = rts[t] + frac * (rts[t - 1] - rts[t]);
                    foundRight = true;
                    break;
                }
            }
            if (!foundRight) return float.NaN;

            return rightRt - leftRt;
        }

        // ══════════════════════════════════════════════════════════════════
        // Feature sub-computations
        // ══════════════════════════════════════════════════════════════════

        private static void AlignToGrid(
            ReadOnlySpan<float> refRts, ReadOnlySpan<float> fragRts,
            ReadOnlySpan<float> fragInts, float[] matrix,
            int fragIdx, int fragCount, int timePoints)
        {
            const float rtTol = 0.01f;
            int fp = 0;
            for (int t = 0; t < timePoints && fp < fragRts.Length; t++)
            {
                float rrt = refRts[t];
                while (fp < fragRts.Length && fragRts[fp] < rrt - rtTol) fp++;
                if (fp < fragRts.Length && MathF.Abs(fragRts[fp] - rrt) <= rtTol)
                    matrix[t * fragCount + fragIdx] = fragInts[fp++];
            }
        }

        private static (float Mean, float Min) PairwiseXicCorrelation(
            float[] matrix, int fc, int timePoints)
        {
            if (fc < 2) return (float.NaN, float.NaN);
            float[] rm = ArrayPool<float>.Shared.Rent(fc);
            float[] rs = ArrayPool<float>.Shared.Rent(fc);
            try
            {
                for (int f = 0; f < fc; f++)
                {
                    double s = 0, s2 = 0;
                    for (int t = 0; t < timePoints; t++) { double v = matrix[t * fc + f]; s += v; s2 += v * v; }
                    double m = s / timePoints; double var = s2 / timePoints - m * m;
                    rm[f] = (float)m; rs[f] = var > 0 ? (float)Math.Sqrt(var) : 0f;
                }
                float sumR = 0f, minR = float.MaxValue; int pairs = 0;
                for (int a = 0; a < fc - 1; a++)
                {
                    if (rs[a] <= 0f) continue;
                    for (int b = a + 1; b < fc; b++)
                    {
                        if (rs[b] <= 0f) continue;
                        double cov = 0;
                        for (int t = 0; t < timePoints; t++)
                            cov += (matrix[t * fc + a] - rm[a]) * (matrix[t * fc + b] - rm[b]);
                        float r = Math.Clamp((float)(cov / timePoints / rs[a] / rs[b]), -1f, 1f);
                        sumR += r; if (r < minR) minR = r; pairs++;
                    }
                }
                return pairs > 0 ? (sumR / pairs, minR) : (float.NaN, float.NaN);
            }
            finally { ArrayPool<float>.Shared.Return(rm); ArrayPool<float>.Shared.Return(rs); }
        }

        /// <summary>
        /// BestFragCorrSum = sum of pairwise correlations for the best fragment.
        /// MedianFragRefCorr = median correlation of all others with the best fragment.
        /// </summary>
        private static (float BestFragCorrSum, float MedianFragRefCorr) ComputeBestFragFeatures(
            float[] matrix, int fc, int timePoints)
        {
            if (fc < 2) return (float.NaN, float.NaN);
            float[] rm = ArrayPool<float>.Shared.Rent(fc);
            float[] rs = ArrayPool<float>.Shared.Rent(fc);
            try
            {
                for (int f = 0; f < fc; f++)
                {
                    double s = 0, s2 = 0;
                    for (int t = 0; t < timePoints; t++) { double v = matrix[t * fc + f]; s += v; s2 += v * v; }
                    double m = s / timePoints; double var = s2 / timePoints - m * m;
                    rm[f] = (float)m; rs[f] = var > 0 ? (float)Math.Sqrt(var) : 0f;
                }
                float[] corrMat = ArrayPool<float>.Shared.Rent(fc * fc);
                try
                {
                    for (int a = 0; a < fc; a++) corrMat[a * fc + a] = 1f;
                    for (int a = 0; a < fc - 1; a++)
                        for (int b = a + 1; b < fc; b++)
                        {
                            float r = 0f;
                            if (rs[a] > 0f && rs[b] > 0f)
                            {
                                double cov = 0;
                                for (int t = 0; t < timePoints; t++)
                                    cov += (matrix[t * fc + a] - rm[a]) * (matrix[t * fc + b] - rm[b]);
                                r = Math.Clamp((float)(cov / timePoints / rs[a] / rs[b]), -1f, 1f);
                            }
                            corrMat[a * fc + b] = r;
                            corrMat[b * fc + a] = r;
                        }
                    int bestIdx = 0; float bestSum = float.MinValue;
                    for (int a = 0; a < fc; a++)
                    {
                        float rowSum = 0f;
                        for (int b = 0; b < fc; b++) if (b != a) rowSum += corrMat[a * fc + b];
                        if (rowSum > bestSum) { bestSum = rowSum; bestIdx = a; }
                    }
                    var refCorrs = new float[fc - 1]; int idx = 0;
                    for (int a = 0; a < fc; a++) if (a != bestIdx) refCorrs[idx++] = corrMat[bestIdx * fc + a];
                    Array.Sort(refCorrs);
                    float median = refCorrs.Length % 2 == 1
                        ? refCorrs[refCorrs.Length / 2]
                        : (refCorrs[refCorrs.Length / 2 - 1] + refCorrs[refCorrs.Length / 2]) * 0.5f;
                    return (bestSum, median);
                }
                finally { ArrayPool<float>.Shared.Return(corrMat); }
            }
            finally { ArrayPool<float>.Shared.Return(rm); ArrayPool<float>.Shared.Return(rs); }
        }

        /// <summary>
        /// Pearson r between MS1 M0 XIC and best fragment XIC. Returns NaN if MS1 unavailable.
        /// Uses the isolation window center as the precursor m/z for MS1 lookup.
        /// </summary>
        private static float ComputeBootstrapMs1Ms2Correlation(
            DiaScanIndex scanIndex,
            int windowId,
            float apexRt,
            ExtractionResult extr,
            int bestFragQueryIdx,
            float ppmTol)
        {
            if (scanIndex.Ms1ScanCount == 0) return float.NaN;
            var bestFragResult = extr.Results[bestFragQueryIdx];
            int ms2Pts = bestFragResult.DataPointCount;
            if (ms2Pts < 3) return float.NaN;

            float rtLo = extr.RtBuffer[bestFragResult.RtBufferOffset];
            float rtHi = extr.RtBuffer[bestFragResult.RtBufferOffset + ms2Pts - 1];

            // Use isolation window center as precursor m/z
            var (wLo, wHi) = scanIndex.GetWindowBounds(windowId);
            float precMz = (wLo + wHi) * 0.5f;
            if (precMz <= 0f) return float.NaN;
            float tol = precMz * ppmTol * 1e-6f;

            var ms1RtList = new List<float>(32);
            var ms1IntList = new List<float>(32);

            for (int si = 0; si < scanIndex.Ms1ScanCount; si++)
            {
                float rt = scanIndex.GetMs1ScanRt(si);
                if (rt < rtLo || rt > rtHi) continue;
                scanIndex.GetMs1ScanPeaks(si, out var mzs, out var ints);
                if (mzs.Length == 0) continue;
                int lo = 0, hi = mzs.Length - 1;
                while (lo < hi) { int mid = (lo + hi) / 2; if (mzs[mid] < precMz - tol) lo = mid + 1; else hi = mid; }
                float bestI = 0f;
                for (int k = lo; k < mzs.Length && mzs[k] <= precMz + tol; k++)
                    if (ints[k] > bestI) bestI = ints[k];
                ms1RtList.Add(rt); ms1IntList.Add(bestI);
            }
            if (ms1RtList.Count < 3) return float.NaN;

            // Align to MS2 RT grid via nearest neighbour
            float[] ms1Aligned = new float[ms2Pts];
            int mp = 0;
            for (int t = 0; t < ms2Pts; t++)
            {
                float rt2 = extr.RtBuffer[bestFragResult.RtBufferOffset + t];
                while (mp < ms1RtList.Count - 1 && Math.Abs(ms1RtList[mp + 1] - rt2) < Math.Abs(ms1RtList[mp] - rt2)) mp++;
                ms1Aligned[t] = ms1IntList[mp];
            }

            ReadOnlySpan<float> ms2Span = new ReadOnlySpan<float>(
                extr.IntensityBuffer, bestFragResult.IntensityBufferOffset, ms2Pts);
            double sm1 = 0, sm2 = 0;
            for (int t = 0; t < ms2Pts; t++) { sm1 += ms1Aligned[t]; sm2 += ms2Span[t]; }
            double m1 = sm1 / ms2Pts, m2 = sm2 / ms2Pts;
            double cov = 0, v1 = 0, v2 = 0;
            for (int t = 0; t < ms2Pts; t++)
            {
                double d1 = ms1Aligned[t] - m1, d2 = ms2Span[t] - m2;
                cov += d1 * d2; v1 += d1 * d1; v2 += d2 * d2;
            }
            if (v1 <= 0 || v2 <= 0) return float.NaN;
            return Math.Clamp((float)(cov / Math.Sqrt(v1 * v2)), -1f, 1f);
        }

        private static float ComputeMeanSignalRatioDev(ReadOnlySpan<float> obs, ReadOnlySpan<float> exp)
        {
            float ot = 0, et = 0;
            for (int f = 0; f < obs.Length; f++) { ot += obs[f]; et += exp[f]; }
            if (ot <= 0 || et <= 0) return float.NaN;
            float s = 0; int n = 0;
            for (int f = 0; f < obs.Length; f++)
            {
                if (obs[f] <= 0 || exp[f] <= 0) continue;
                s += MathF.Abs(MathF.Log2((obs[f] / ot) / (exp[f] / et))); n++;
            }
            return n >= 3 ? s / n : float.NaN;
        }

        private static float ComputeMassErrorPpm(
            SliceGroup g, FragmentQuery[] queries, ExtractionResult extr,
            DiaScanIndex scanIndex, float apexRt, float ppmTol)
        {
            if (!scanIndex.TryGetScanRangeForWindow(g.WindowId, out int wS, out int wC)) return ppmTol;
            int ai = -1;
            for (int si = wS; si < wS + wC; si++)
                if (MathF.Abs(scanIndex.GetScanRt(si) - apexRt) < 1e-5f) { ai = si; break; }
            if (ai < 0) return ppmTol;
            ReadOnlySpan<float> mzs = scanIndex.GetScanMzSpan(ai);
            if (mzs.Length == 0) return ppmTol;
            double eS = 0; int eN = 0;
            for (int qi = g.QueryOffset; qi < g.QueryOffset + g.QueryCount; qi++)
            {
                var res = extr.Results[qi];
                if (res.DataPointCount == 0) continue;
                bool hit = false;
                for (int k = 0; k < res.DataPointCount; k++)
                    if (extr.RtBuffer[res.RtBufferOffset + k] == apexRt &&
                       extr.IntensityBuffer[res.IntensityBufferOffset + k] > 0) { hit = true; break; }
                if (!hit) continue;
                float tgt = queries[qi].TargetMz;
                int lo = 0, hi = mzs.Length - 1;
                while (lo < hi) { int mid = (lo + hi) / 2; if (mzs[mid] < tgt) lo = mid + 1; else hi = mid; }
                float best = mzs[lo];
                if (lo > 0 && MathF.Abs(mzs[lo - 1] - tgt) < MathF.Abs(best - tgt)) best = mzs[lo - 1];
                float ppm = MathF.Abs((best - tgt) / tgt * 1e6f);
                if (ppm <= ppmTol) { eS += ppm; eN++; }
            }
            return eN > 0 ? (float)(eS / eN) : ppmTol;
        }

        private static float MedianOfSpan(Span<float> src, int n)
        {
            if (n == 0) return 0f;
            float[] r = ArrayPool<float>.Shared.Rent(n);
            try
            {
                src.Slice(0, n).CopyTo(r); Array.Sort(r, 0, n);
                return n % 2 == 0 ? (r[n / 2 - 1] + r[n / 2]) * 0.5f : r[n / 2];
            }
            finally { ArrayPool<float>.Shared.Return(r); }
        }

        private static float MedianList(List<float> vals)
        {
            if (vals.Count == 0) return 0f;
            var sorted = vals.ToArray(); Array.Sort(sorted);
            int n = sorted.Length;
            return n % 2 == 0 ? (sorted[n / 2 - 1] + sorted[n / 2]) * 0.5f : sorted[n / 2];
        }

        // ══════════════════════════════════════════════════════════════════
        // Iterative Fisher LDA
        // ══════════════════════════════════════════════════════════════════

        private static BootstrapLdaModel FitIterativeLda(
            float[][] features, bool[] isDecoy,
            float[]? warmStartWeights, Action<string> log)
        {
            int n = features.Length;
            bool[] labels = (bool[])isDecoy.Clone();

            if (warmStartWeights != null && warmStartWeights.Length == FeatureCount)
            {
                float[] s = new float[n];
                for (int i = 0; i < n; i++) for (int j = 0; j < FeatureCount; j++) s[i] += features[i][j] * warmStartWeights[j];
                double tS = 0, dS = 0; int tN = 0, dN = 0;
                for (int i = 0; i < n; i++) if (!isDecoy[i]) { tS += s[i]; tN++; } else { dS += s[i]; dN++; }
                float wt = tN > 0 && dN > 0 ? (float)((tS / tN + dS / dN) * 0.5) : 0f;
                for (int i = 0; i < n; i++) labels[i] = s[i] < wt;
                log("Using warm-start weights for initial labelling.");
            }

            float[] weights = new float[FeatureCount];
            float threshold = 0f, prevBound = float.MaxValue;
            int iter = 0;

            for (iter = 0; iter < MaxIterations; iter++)
            {
                FitFisherLda(features, labels, out weights, out threshold);
                float[] scores = new float[n];
                for (int i = 0; i < n; i++) scores[i] = Dot(features[i], weights);

                (float tM, float tS) = MeanSigma(scores, labels, true);
                (float dM, _) = MeanSigma(scores, labels, false);
                float boundary = (tM + dM) * 0.5f;
                float movement = MathF.Abs(boundary - prevBound);
                prevBound = boundary; threshold = boundary;

                int nT = 0, nD = 0;
                for (int i = 0; i < n; i++) { if (!labels[i]) nT++; else nD++; }
                log($"Iter {iter + 1,2}: boundary={boundary:F4}  tMean={tM:F4}  dMean={dM:F4}  " +
                    $"tSigma={tS:F4}  sep={(tM - dM) / Math.Max(tS, 1e-6f):F2}  nT={nT}  nD={nD}  δ={movement:F5}");

                if (movement < ConvergenceEpsilon && iter > 0)
                { log($"Converged after {iter + 1} iterations (δ={movement:F5})."); break; }

                for (int i = 0; i < n; i++) labels[i] = scores[i] < boundary;
            }

            float[] fs = new float[n];
            for (int i = 0; i < n; i++) fs[i] = Dot(features[i], weights);
            (float ftM, float ftS) = MeanSigma(fs, labels, true);
            (float fdM, _) = MeanSigma(fs, labels, false);

            log($"Final: threshold={threshold:F4}  tMean={ftM:F4}  dMean={fdM:F4}  " +
                $"tSigma={ftS:F4}  sep={(ftM - fdM) / Math.Max(ftS, 1e-6f):F2}");
            for (int j = 0; j < FeatureCount; j++)
                log($"  w[{j}] {FeatureNames[j],-22} = {weights[j]:+0.0000;-0.0000}");

            return new BootstrapLdaModel(weights, threshold, ftM, ftS, fdM, iter + 1);
        }

        internal static void FitFisherLda(
            float[][] features, bool[] labels,
            out float[] weights, out float threshold)
        {
            int n = features.Length, d = FeatureCount;
            double[] mu0 = new double[d], mu1 = new double[d];
            double[] v0 = new double[d], v1 = new double[d];
            int n0 = 0, n1 = 0;
            for (int i = 0; i < n; i++) { if (labels[i]) n1++; else n0++; }
            for (int i = 0; i < n; i++) { var mu = labels[i] ? mu1 : mu0; for (int j = 0; j < d; j++) mu[j] += features[i][j]; }
            if (n0 > 0) for (int j = 0; j < d; j++) mu0[j] /= n0;
            if (n1 > 0) for (int j = 0; j < d; j++) mu1[j] /= n1;
            for (int i = 0; i < n; i++)
            {
                var mu = labels[i] ? mu1 : mu0; var v = labels[i] ? v1 : v0;
                for (int j = 0; j < d; j++) { double diff = features[i][j] - mu[j]; v[j] += diff * diff; }
            }
            weights = new float[d]; double l1 = 0;
            for (int j = 0; j < d; j++)
            {
                double sw = Math.Max((v0[j] + v1[j]) / Math.Max(n - 2, 1), 1e-8);
                weights[j] = (float)((mu0[j] - mu1[j]) / sw); l1 += Math.Abs(weights[j]);
            }
            if (l1 > 0) for (int j = 0; j < d; j++) weights[j] = (float)(weights[j] / l1);
            double p0 = 0, p1 = 0;
            for (int j = 0; j < d; j++) { p0 += weights[j] * mu0[j]; p1 += weights[j] * mu1[j]; }
            threshold = (float)((p0 + p1) * 0.5);
        }

        private static (float Mean, float Sigma) MeanSigma(float[] s, bool[] l, bool isTarget)
        {
            double sum = 0; int cnt = 0;
            for (int i = 0; i < s.Length; i++) if (l[i] != isTarget) { sum += s[i]; cnt++; }
            if (cnt == 0) return (0f, 1f);
            double mean = sum / cnt, var = 0;
            for (int i = 0; i < s.Length; i++) if (l[i] != isTarget) { double d = s[i] - mean; var += d * d; }
            return ((float)mean, (float)Math.Sqrt(var / Math.Max(cnt - 1, 1)));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float Dot(float[] f, float[] w)
        { float s = 0f; for (int j = 0; j < FeatureCount; j++) s += f[j] * w[j]; return s; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static float NdpScore(ReadOnlySpan<float> obs, ReadOnlySpan<float> exp)
        {
            float dot = 0f, no = 0f, ne = 0f;
            int len = Math.Min(obs.Length, exp.Length);
            for (int i = 0; i < len; i++) { dot += obs[i] * exp[i]; no += obs[i] * obs[i]; ne += exp[i] * exp[i]; }
            if (no <= 0f || ne <= 0f) return 0f;
            return Math.Clamp(dot / (MathF.Sqrt(no) * MathF.Sqrt(ne)), 0f, 1f);
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // BootstrapCalibrationResult
    // ══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Result of a full bidirectional bootstrap calibration run.
    /// </summary>
    public sealed class BootstrapCalibrationResult
    {
        /// <summary>
        /// Non-linear LOWESS calibration model fitted on all pooled anchors.
        /// Primary calibration artifact. Null if fewer than 50 anchors were collected.
        /// </summary>
        public LowessRtModel? LowessModel { get; }

        /// <summary>Global OLS line fitted on all pooled anchors (diagnostics / fallback).</summary>
        public BootstrapRtLine? GlobalOlsLine { get; }

        /// <summary>All anchors pooled from seed + left arm + right arm.</summary>
        public IReadOnlyList<BootstrapAnchor> AllAnchors { get; }

        /// <summary>Seed step result. Null only if the seed itself failed.</summary>
        public BootstrapSliceResult? Seed { get; }

        /// <summary>Left-arm step results in order (earliest step first).</summary>
        public IReadOnlyList<BootstrapSliceResult> LeftSteps { get; }

        /// <summary>Right-arm step results in order (earliest step first).</summary>
        public IReadOnlyList<BootstrapSliceResult> RightSteps { get; }

        internal BootstrapCalibrationResult(
            LowessRtModel? lowess,
            BootstrapRtLine? globalOls,
            List<BootstrapAnchor> allAnchors,
            BootstrapSliceResult? seed,
            List<BootstrapSliceResult> leftSteps,
            List<BootstrapSliceResult> rightSteps)
        {
            LowessModel = lowess;
            GlobalOlsLine = globalOls;
            AllAnchors = allAnchors;
            Seed = seed;
            LeftSteps = leftSteps;
            RightSteps = rightSteps;
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // BootstrapSliceResult
    // ══════════════════════════════════════════════════════════════════════

    public sealed class BootstrapSliceResult
    {
        public BootstrapLdaModel Model { get; }
        public BootstrapRtLine? RtLine { get; }
        public IReadOnlyList<BootstrapAnchor> Anchors { get; }
        public float SliceLo { get; }
        public float SliceHi { get; }
        /// <summary>Median FWHM of confident target peaks in this slice (minutes). 0 if not measured.</summary>
        public float MedianPeakFwhmMinutes { get; }

        internal BootstrapSliceResult(
            BootstrapLdaModel m, BootstrapRtLine? l,
            List<BootstrapAnchor> a, float lo, float hi, float fwhm)
        { Model = m; RtLine = l; Anchors = a; SliceLo = lo; SliceHi = hi; MedianPeakFwhmMinutes = fwhm; }

        /// <summary>Shifts the slice left by stepPct of total RT range. Returns null at left edge.</summary>
        public (float Lo, float Hi)? ShiftLeft(DiaScanIndex scanIndex, float stepPct = DiaBootstrapCalibrator.StepPct)
        {
            int n = scanIndex.ScanCount; if (n == 0) return null;
            float[] buf = ArrayPool<float>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) buf[i] = scanIndex.GetScanRt(i);
                Array.Sort(buf, 0, n);
                float shift = (buf[n - 1] - buf[0]) * stepPct;
                float newLo = SliceLo - shift, newHi = SliceHi - shift;
                if (newHi <= buf[0]) return null;
                return (Math.Max(newLo, buf[0]), newHi);
            }
            finally { ArrayPool<float>.Shared.Return(buf); }
        }

        /// <summary>Shifts the slice right by stepPct of total RT range. Returns null at right edge.</summary>
        public (float Lo, float Hi)? ShiftRight(DiaScanIndex scanIndex, float stepPct = DiaBootstrapCalibrator.StepPct)
        {
            int n = scanIndex.ScanCount; if (n == 0) return null;
            float[] buf = ArrayPool<float>.Shared.Rent(n);
            try
            {
                for (int i = 0; i < n; i++) buf[i] = scanIndex.GetScanRt(i);
                Array.Sort(buf, 0, n);
                float shift = (buf[n - 1] - buf[0]) * stepPct;
                float newLo = SliceLo + shift, newHi = SliceHi + shift;
                if (newLo >= buf[n - 1]) return null;
                return (newLo, Math.Min(newHi, buf[n - 1]));
            }
            finally { ArrayPool<float>.Shared.Return(buf); }
        }
    }

    // ══════════════════════════════════════════════════════════════════════
    // BootstrapAnchor
    // ══════════════════════════════════════════════════════════════════════

    public readonly struct BootstrapAnchor
    {
        public readonly float LibraryRt;
        public readonly float ObservedApexRt;
        public readonly float LdaScore;
        public BootstrapAnchor(float lib, float obs, float score)
        { LibraryRt = lib; ObservedApexRt = obs; LdaScore = score; }
    }

    // ══════════════════════════════════════════════════════════════════════
    // BootstrapRtLine
    // ══════════════════════════════════════════════════════════════════════

    public sealed class BootstrapRtLine
    {
        public float Slope { get; }
        public float Intercept { get; }
        public float RSquared { get; }
        public float ResidualSigma { get; }
        /// <summary>HalfWidth = max(ResidualSigma × 2, fwhmFloor). Never goes below fwhmFloor.</summary>
        public float HalfWidth { get; }
        public int AnchorCount { get; }

        private BootstrapRtLine(float s, float i, float r2, float sig, float hw, int n)
        { Slope = s; Intercept = i; RSquared = r2; ResidualSigma = sig; HalfWidth = hw; AnchorCount = n; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public float PredictObservedRt(float libraryRt) => Slope * libraryRt + Intercept;

        /// <summary>OLS fit. HalfWidth = max(σ×2, fwhmFloor).</summary>
        public static BootstrapRtLine? Fit(
            IReadOnlyList<BootstrapAnchor> anchors, float fwhmFloor = 0f)
        {
            double sX = 0, sY = 0, sXY = 0, sX2 = 0; int n = 0;
            foreach (var a in anchors)
            {
                if (a.LibraryRt <= 0f || a.ObservedApexRt <= 0f) continue;
                sX += a.LibraryRt; sY += a.ObservedApexRt;
                sXY += a.LibraryRt * a.ObservedApexRt;
                sX2 += a.LibraryRt * a.LibraryRt; n++;
            }
            if (n < 5) return null;
            double den = n * sX2 - sX * sX;
            if (Math.Abs(den) < 1e-12) return null;
            float slope = (float)((n * sXY - sX * sY) / den);
            float intercept = (float)((sY - slope * sX) / n);
            return ComputeStats(anchors, slope, intercept, n, fwhmFloor);
        }

        /// <summary>
        /// Fits only the intercept with a fixed slope (slope rescue).
        /// intercept = mean(observedRt) - slope × mean(libraryRt)
        /// </summary>
        public static BootstrapRtLine? FitWithFixedSlope(
            IReadOnlyList<BootstrapAnchor> anchors, float slope, float fwhmFloor = 0f)
        {
            double sX = 0, sY = 0; int n = 0;
            foreach (var a in anchors)
            {
                if (a.LibraryRt <= 0f || a.ObservedApexRt <= 0f) continue;
                sX += a.LibraryRt; sY += a.ObservedApexRt; n++;
            }
            if (n < 3) return null;
            float intercept = (float)(sY / n - slope * (sX / n));
            return ComputeStats(anchors, slope, intercept, n, fwhmFloor);
        }

        private static BootstrapRtLine ComputeStats(
            IReadOnlyList<BootstrapAnchor> anchors,
            float slope, float intercept, int n, float fwhmFloor)
        {
            double meanY = 0, ssTot = 0, ssRes = 0;
            foreach (var a in anchors)
                if (a.LibraryRt > 0f && a.ObservedApexRt > 0f) meanY += a.ObservedApexRt;
            meanY /= n;
            foreach (var a in anchors)
            {
                if (a.LibraryRt <= 0f || a.ObservedApexRt <= 0f) continue;
                double pred = slope * a.LibraryRt + intercept;
                double res = a.ObservedApexRt - pred;
                ssRes += res * res;
                double dev = a.ObservedApexRt - meanY;
                ssTot += dev * dev;
            }
            float r2 = ssTot > 0 ? (float)(1.0 - ssRes / ssTot) : 0f;
            float sigma = n > 2 ? (float)Math.Sqrt(ssRes / (n - 2)) : 0f;
            float hw = Math.Max(sigma * 2f, fwhmFloor);
            return new BootstrapRtLine(slope, intercept, r2, sigma, hw, n);
        }

        public override string ToString() =>
            $"obs = {Slope:F4}×lib + {Intercept:F4}  " +
            $"R²={RSquared:F4}  σ={ResidualSigma:F4} min  hw={HalfWidth:F4} min  n={AnchorCount}";
    }

    // ══════════════════════════════════════════════════════════════════════
    // BootstrapLdaModel
    // ══════════════════════════════════════════════════════════════════════

    public sealed class BootstrapLdaModel
    {
        public float[] Weights { get; }
        public float Threshold { get; }
        public float TargetMean { get; }
        public float TargetSigma { get; }
        public float DecoyMean { get; }
        public int IterationsRun { get; }

        public float SeparationRatio =>
            TargetSigma > 0f ? (TargetMean - DecoyMean) / TargetSigma : 0f;

        internal BootstrapLdaModel(float[] w, float thr, float tM, float tS, float dM, int it)
        { Weights = w; Threshold = thr; TargetMean = tM; TargetSigma = tS; DecoyMean = dM; IterationsRun = it; }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public float Score(ReadOnlySpan<float> fv)
        {
            float s = 0f; int len = Math.Min(fv.Length, Weights.Length);
            for (int i = 0; i < len; i++) s += fv[i] * Weights[i]; return s;
        }

        public bool IsTarget(float score) => score >= Threshold;

        public override string ToString() =>
            $"LDA  threshold={Threshold:F4}  tMean={TargetMean:F4}  " +
            $"dMean={DecoyMean:F4}  tSigma={TargetSigma:F4}  sep={SeparationRatio:F2}  iters={IterationsRun}";
    }

    // ── Internal ──────────────────────────────────────────────────────────

    internal readonly struct SliceGroup
    {
        public readonly int PrecursorIndex, QueryOffset, QueryCount, WindowId;
        public SliceGroup(int pi, int qo, int qc, int wid)
        { PrecursorIndex = pi; QueryOffset = qo; QueryCount = qc; WindowId = wid; }
    }
}