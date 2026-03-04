// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Phase 15, Prompt 3: Non-linear calibration models
// Placement: MassSpectrometry/Dia/Calibration/PiecewiseLinearRtModel.cs

using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Dia.Calibration
{
    /// <summary>
    /// Enum identifying which calibration model type is in use.
    /// Added in Phase 15 to support automatic model selection during iterative calibration.
    /// </summary>
    public enum RtCalibrationModelType
    {
        /// <summary>Current production model: y = slope * x + intercept.</summary>
        Linear,

        /// <summary>Monotonic piecewise-linear with N knots for gradient-shape mismatches.</summary>
        PiecewiseLinear,

        /// <summary>Locally-weighted scatterplot smoothing with isotonic enforcement.</summary>
        Lowess
    }

    /// <summary>
    /// Internal interface for RT calibration models used during iterative calibration.
    /// 
    /// Because <see cref="RtCalibrationModel"/> is sealed and immutable, non-linear models
    /// cannot inherit from it. Instead, they implement this interface for use within the
    /// iterative calibration loop, and expose a <see cref="ToRtCalibrationModel"/> method
    /// to produce a standard <see cref="RtCalibrationModel"/> for downstream consumers
    /// (query generation, FDR, etc.) that expect the production type.
    /// </summary>
    public interface IRtCalibrationModel
    {
        /// <summary>Predicts observed RT (minutes) from a library RT / iRT value.</summary>
        double ToMinutes(double libraryRtOrIrt);

        /// <summary>Global residual σ in minutes.</summary>
        double SigmaMinutes { get; }

        /// <summary>R² goodness-of-fit.</summary>
        double RSquared { get; }

        /// <summary>Number of anchors used in fitting.</summary>
        int AnchorCount { get; }

        /// <summary>Whether this model is reliable enough for use.</summary>
        bool IsReliable { get; }

        /// <summary>Model type identifier.</summary>
        RtCalibrationModelType ModelType { get; }

        /// <summary>Global slope (for logging / backward compat).</summary>
        double Slope { get; }

        /// <summary>Global intercept (for logging / backward compat).</summary>
        double Intercept { get; }

        /// <summary>
        /// Returns the local σ for a given library RT. For linear models this is
        /// just the global σ; for piecewise/LOWESS models it varies by RT region.
        /// </summary>
        double GetLocalSigma(double libraryRtOrIrt);

        /// <summary>
        /// Produces a standard <see cref="RtCalibrationModel"/> for downstream code
        /// that expects the sealed production type.
        /// </summary>
        RtCalibrationModel ToRtCalibrationModel();
    }

    /// <summary>
    /// Wraps an existing <see cref="RtCalibrationModel"/> to implement <see cref="IRtCalibrationModel"/>.
    /// Zero overhead — just delegates all calls.
    /// </summary>
    public sealed class LinearRtModelWrapper : IRtCalibrationModel
    {
        private readonly RtCalibrationModel _model;

        public LinearRtModelWrapper(RtCalibrationModel model)
        {
            _model = model ?? throw new ArgumentNullException(nameof(model));
        }

        public double ToMinutes(double libraryRtOrIrt) => _model.ToMinutes(libraryRtOrIrt);
        public double SigmaMinutes => _model.SigmaMinutes;
        public double RSquared => _model.RSquared;
        public int AnchorCount => _model.AnchorCount;
        public bool IsReliable => _model.IsReliable;
        public RtCalibrationModelType ModelType => RtCalibrationModelType.Linear;
        public double Slope => _model.Slope;
        public double Intercept => _model.Intercept;
        public double GetLocalSigma(double libraryRtOrIrt) => _model.SigmaMinutes;
        public RtCalibrationModel ToRtCalibrationModel() => _model;
    }

    /// <summary>
    /// Piecewise-linear RT calibration model that divides the library RT range into N segments,
    /// fits a separate linear model per segment, and enforces monotonicity at segment boundaries.
    ///
    /// This handles gradient-shape mismatches between library and observed RTs that a single
    /// global linear model cannot capture — e.g., libraries from different instruments, predicted
    /// RTs from deep-learning models, or localized RT shifts from column degradation.
    ///
    /// Because <see cref="RtCalibrationModel"/> is sealed, this class implements
    /// <see cref="IRtCalibrationModel"/> instead.
    /// </summary>
    public sealed class PiecewiseLinearRtModel : IRtCalibrationModel
    {
        // ── Segment boundary data ──────────────────────────────────────────

        /// <summary>Library RT values at each knot (segment boundary), sorted ascending. Length = NumSegments + 1.</summary>
        public float[] KnotLibraryRts { get; }

        /// <summary>Predicted (observed) RT at each knot after monotonicity enforcement. Length = NumSegments + 1.</summary>
        public float[] KnotPredictedRts { get; }

        /// <summary>Per-segment slopes. Length = NumSegments.</summary>
        public float[] SegmentSlopes { get; }

        /// <summary>Per-segment intercepts. Length = NumSegments.</summary>
        public float[] SegmentIntercepts { get; }

        /// <summary>Per-segment residual σ in minutes. Length = NumSegments.</summary>
        public float[] SegmentSigmas { get; }

        /// <summary>Number of anchors used in each segment. Length = NumSegments.</summary>
        public int[] SegmentAnchorCounts { get; }

        /// <summary>Number of segments after merging for monotonicity.</summary>
        public int NumSegments => SegmentSlopes.Length;

        // ── IRtCalibrationModel ────────────────────────────────────────────

        public double SigmaMinutes { get; }
        public double RSquared { get; }
        public int AnchorCount { get; }
        public bool IsReliable => RSquared >= RtCalibrationModel.MinReliableRSquared
                               && AnchorCount >= RtCalibrationModel.MinReliableAnchors;
        public RtCalibrationModelType ModelType => RtCalibrationModelType.PiecewiseLinear;
        public double Slope { get; }
        public double Intercept { get; }

        // ── Construction ───────────────────────────────────────────────────

        private PiecewiseLinearRtModel(
            float[] knotLibraryRts,
            float[] knotPredictedRts,
            float[] segmentSlopes,
            float[] segmentIntercepts,
            float[] segmentSigmas,
            int[] segmentAnchorCounts,
            double globalSlope,
            double globalIntercept,
            double globalSigma,
            double globalRSquared,
            int totalAnchorCount)
        {
            KnotLibraryRts = knotLibraryRts;
            KnotPredictedRts = knotPredictedRts;
            SegmentSlopes = segmentSlopes;
            SegmentIntercepts = segmentIntercepts;
            SegmentSigmas = segmentSigmas;
            SegmentAnchorCounts = segmentAnchorCounts;
            Slope = globalSlope;
            Intercept = globalIntercept;
            SigmaMinutes = globalSigma;
            RSquared = globalRSquared;
            AnchorCount = totalAnchorCount;
        }

        // ── Prediction ─────────────────────────────────────────────────────

        public double ToMinutes(double libraryRtOrIrt)
        {
            if (libraryRtOrIrt <= KnotLibraryRts[0])
                return SegmentSlopes[0] * libraryRtOrIrt + SegmentIntercepts[0];

            if (libraryRtOrIrt >= KnotLibraryRts[KnotLibraryRts.Length - 1])
            {
                int last = NumSegments - 1;
                return SegmentSlopes[last] * libraryRtOrIrt + SegmentIntercepts[last];
            }

            int segIdx = FindSegment(libraryRtOrIrt);
            float leftRt = KnotLibraryRts[segIdx];
            float rightRt = KnotLibraryRts[segIdx + 1];
            float leftPred = KnotPredictedRts[segIdx];
            float rightPred = KnotPredictedRts[segIdx + 1];

            double fraction = (libraryRtOrIrt - leftRt) / (rightRt - leftRt);
            return leftPred + fraction * (rightPred - leftPred);
        }

        public double GetLocalSigma(double libraryRtOrIrt)
        {
            if (libraryRtOrIrt <= KnotLibraryRts[0])
                return SegmentSigmas[0];
            if (libraryRtOrIrt >= KnotLibraryRts[KnotLibraryRts.Length - 1])
                return SegmentSigmas[NumSegments - 1];
            return SegmentSigmas[FindSegment(libraryRtOrIrt)];
        }

        public RtCalibrationModel ToRtCalibrationModel()
        {
            double slope = Math.Abs(Slope) < 1e-12 ? 1.0 : Slope;
            return new RtCalibrationModel(
                slope: slope,
                intercept: Intercept,
                sigmaMinutes: SigmaMinutes,
                rSquared: RSquared,
                anchorCount: AnchorCount);
        }

        // ── Fitting ────────────────────────────────────────────────────────

        /// <summary>
        /// Fits a piecewise-linear model. Returns null if fewer than 50 anchors.
        /// </summary>
        public static PiecewiseLinearRtModel Fit(
            double[] libraryRts,
            double[] observedRts,
            int numSegments = 8)
        {
            int n = libraryRts.Length;
            if (n < 50) return null;

            var indices = Enumerable.Range(0, n).ToArray();
            var libSorted = new double[n];
            var obsSorted = new double[n];
            Array.Copy(libraryRts, libSorted, n);
            Array.Sort(libSorted, indices);
            for (int i = 0; i < n; i++)
                obsSorted[i] = observedRts[indices[i]];

            int maxSegments = Math.Max(2, n / 10);
            numSegments = Math.Min(numSegments, maxSegments);

            var segments = SplitIntoSegments(libSorted, obsSorted, n, numSegments);
            FitSegmentOLS(segments);
            segments = EnforceMonotonicity(segments);
            return BuildModel(segments, libSorted, obsSorted, n);
        }

        // ── Private helpers ────────────────────────────────────────────────

        private int FindSegment(double libraryRt)
        {
            int lo = 0, hi = KnotLibraryRts.Length - 2;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (libraryRt < KnotLibraryRts[mid + 1]) hi = mid;
                else lo = mid + 1;
            }
            return lo;
        }

        private class SegmentData
        {
            public double LibRtMin, LibRtMax;
            public List<double> LibRts, ObsRts;
            public double Slope, Intercept, Sigma;
        }

        private static List<SegmentData> SplitIntoSegments(
            double[] sortedLib, double[] sortedObs, int n, int numSeg)
        {
            var segments = new List<SegmentData>(numSeg);
            int baseSize = n / numSeg;
            int remainder = n % numSeg;
            int offset = 0;

            for (int s = 0; s < numSeg; s++)
            {
                int count = baseSize + (s < remainder ? 1 : 0);
                var seg = new SegmentData
                {
                    LibRts = new List<double>(count),
                    ObsRts = new List<double>(count)
                };
                for (int i = offset; i < offset + count; i++)
                {
                    seg.LibRts.Add(sortedLib[i]);
                    seg.ObsRts.Add(sortedObs[i]);
                }
                seg.LibRtMin = seg.LibRts[0];
                seg.LibRtMax = seg.LibRts[seg.LibRts.Count - 1];
                segments.Add(seg);
                offset += count;
            }
            return segments;
        }

        private static void FitSegmentOLS(List<SegmentData> segments)
        {
            foreach (var seg in segments)
            {
                if (seg.LibRts.Count < 2)
                {
                    seg.Slope = 1.0;
                    seg.Intercept = seg.ObsRts[0] - seg.LibRts[0];
                }
                else
                {
                    FitOLS(seg.LibRts, seg.ObsRts, out double s, out double ic);
                    seg.Slope = s;
                    seg.Intercept = ic;
                }

                double sumSq = 0;
                for (int i = 0; i < seg.LibRts.Count; i++)
                {
                    double resid = seg.ObsRts[i] - (seg.Slope * seg.LibRts[i] + seg.Intercept);
                    sumSq += resid * resid;
                }
                seg.Sigma = Math.Sqrt(sumSq / seg.LibRts.Count);
            }
        }

        private static List<SegmentData> EnforceMonotonicity(List<SegmentData> segments)
        {
            bool merged = true;
            while (merged)
            {
                merged = false;
                for (int i = 0; i < segments.Count - 1; i++)
                {
                    double predRight = segments[i].Slope * segments[i].LibRtMax + segments[i].Intercept;
                    double predLeftNext = segments[i + 1].Slope * segments[i + 1].LibRtMin + segments[i + 1].Intercept;

                    if (predRight > predLeftNext)
                    {
                        var combined = new SegmentData
                        {
                            LibRts = new List<double>(segments[i].LibRts.Count + segments[i + 1].LibRts.Count),
                            ObsRts = new List<double>(segments[i].ObsRts.Count + segments[i + 1].ObsRts.Count)
                        };
                        combined.LibRts.AddRange(segments[i].LibRts);
                        combined.LibRts.AddRange(segments[i + 1].LibRts);
                        combined.ObsRts.AddRange(segments[i].ObsRts);
                        combined.ObsRts.AddRange(segments[i + 1].ObsRts);
                        combined.LibRtMin = combined.LibRts[0];
                        combined.LibRtMax = combined.LibRts[combined.LibRts.Count - 1];

                        FitOLS(combined.LibRts, combined.ObsRts, out double s, out double ic);
                        combined.Slope = s;
                        combined.Intercept = ic;
                        double sumSq = 0;
                        for (int j = 0; j < combined.LibRts.Count; j++)
                        {
                            double r = combined.ObsRts[j] - (s * combined.LibRts[j] + ic);
                            sumSq += r * r;
                        }
                        combined.Sigma = Math.Sqrt(sumSq / combined.LibRts.Count);

                        segments[i] = combined;
                        segments.RemoveAt(i + 1);
                        merged = true;
                        break;
                    }
                }
            }
            return segments;
        }

        private static PiecewiseLinearRtModel BuildModel(
            List<SegmentData> segments, double[] sortedLib, double[] sortedObs, int n)
        {
            int numSeg = segments.Count;
            var knotLib = new float[numSeg + 1];
            var knotPred = new float[numSeg + 1];
            var slopes = new float[numSeg];
            var intercepts = new float[numSeg];
            var sigmas = new float[numSeg];
            var counts = new int[numSeg];

            knotLib[0] = (float)segments[0].LibRtMin;
            knotPred[0] = (float)(segments[0].Slope * segments[0].LibRtMin + segments[0].Intercept);

            for (int s = 0; s < numSeg; s++)
            {
                slopes[s] = (float)segments[s].Slope;
                intercepts[s] = (float)segments[s].Intercept;
                sigmas[s] = (float)segments[s].Sigma;
                counts[s] = segments[s].LibRts.Count;
                knotLib[s + 1] = (float)segments[s].LibRtMax;
                knotPred[s + 1] = (float)(segments[s].Slope * segments[s].LibRtMax + segments[s].Intercept);
            }

            for (int i = 1; i < knotPred.Length; i++)
                if (knotPred[i] < knotPred[i - 1]) knotPred[i] = knotPred[i - 1];

            FitOLS(sortedLib, sortedObs, n, out double gSlope, out double gIntercept);

            // Temporary model for computing piecewise predictions → global σ, R²
            var temp = new PiecewiseLinearRtModel(
                knotLib, knotPred, slopes, intercepts, sigmas, counts,
                gSlope, gIntercept, 0, 0, n);

            double ssRes = 0, ssTot = 0, meanObs = 0;
            for (int i = 0; i < n; i++) meanObs += sortedObs[i];
            meanObs /= n;
            for (int i = 0; i < n; i++)
            {
                double pred = temp.ToMinutes(sortedLib[i]);
                double resid = sortedObs[i] - pred;
                ssRes += resid * resid;
                double dev = sortedObs[i] - meanObs;
                ssTot += dev * dev;
            }

            double gSigma = Math.Sqrt(ssRes / n);
            double gRSq = ssTot > 0 ? 1.0 - ssRes / ssTot : 0;

            return new PiecewiseLinearRtModel(
                knotLib, knotPred, slopes, intercepts, sigmas, counts,
                gSlope, gIntercept, gSigma, gRSq, n);
        }

        private static void FitOLS(IList<double> xs, IList<double> ys, out double slope, out double intercept)
        {
            int n = xs.Count;
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
            double d = n * sxx - sx * sx;
            if (Math.Abs(d) < 1e-12) { slope = 1.0; intercept = (sy - sx) / n; return; }
            slope = (n * sxy - sx * sy) / d;
            intercept = (sy - slope * sx) / n;
        }

        private static void FitOLS(double[] xs, double[] ys, int n, out double slope, out double intercept)
        {
            double sx = 0, sy = 0, sxx = 0, sxy = 0;
            for (int i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
            double d = n * sxx - sx * sx;
            if (Math.Abs(d) < 1e-12) { slope = 1.0; intercept = (sy - sx) / n; return; }
            slope = (n * sxy - sx * sy) / d;
            intercept = (sy - slope * sx) / n;
        }
    }
}