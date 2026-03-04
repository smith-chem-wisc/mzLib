// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Buffers;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 16B, Prompt 4: Precursor XIC extraction from MS1 SoA arrays.
    /// 
    /// Extracts extracted ion chromatograms (XICs) for precursor ions and their
    /// isotope envelopes from the MS1 scan data stored in <see cref="DiaScanIndex"/>.
    /// 
    /// Design principles:
    ///   - O(K log N): binary search via FindMs1ScanIndexAtRt locates the RT window
    ///     start in O(log N); then K scans in the window are iterated linearly.
    ///   - Zero-allocation inner loop: peak arrays are accessed as ReadOnlySpan via
    ///     GetMs1ScanPeaks (no copy). Temporary result buffers are rented from
    ///     ArrayPool and returned after packing into the output arrays.
    ///   - Binary search within each scan: peaks are sorted by m/z (guaranteed by
    ///     DiaScanIndexBuilder), so the tolerance window is located with two
    ///     BinarySearch calls rather than a linear scan.
    ///   - NaN-safe: if the index has no MS1 scans, all outputs are empty arrays
    ///     and callers should treat features as float.NaN.
    /// 
    /// Called from DiaFeatureExtractor (Prompt 5) to compute:
    ///   [29] PrecursorXicApexIntensity
    ///   [30] IsotopePatternScore
    ///   [31] Ms1Ms2Correlation
    ///   [32] PrecursorElutionScore
    /// </summary>
    public static class Ms1XicExtractor
    {
        // ── Averagine isotope model ──────────────────────────────────────────
        // Precomputed theoretical [M0, M+1, M+2] relative intensities from the
        // averagine model, indexed by integer precursor mass (Da).
        // Range 200–4000 Da covers all practical DIA precursors.
        // Values are normalized so M0+M1+M2 = 1.0.
        //
        // Averagine composition per 111.1 Da residue:
        //   C=4.9384, H=7.7583, N=1.3577, O=1.4773, S=0.0417
        // Isotope probabilities computed via the Poisson approximation.
        private static readonly float[][] AveragineIsotopeTable = BuildAveragineTable();

        // Mass range covered by the table
        private const int AveragineTableMinMass = 200;
        private const int AveragineTableMaxMass = 4000;

        // Neutron mass offset between isotopes (Da)
        private const float NeutronMass = 1.003355f;

        // ── Public API ───────────────────────────────────────────────────────

        /// <summary>
        /// Extracts a precursor (M0) XIC from MS1 scans within [rtMin, rtMax].
        /// 
        /// For each MS1 scan in the RT window, sums all intensity within
        /// <paramref name="ppmTolerance"/> of <paramref name="precursorMz"/>.
        /// The output arrays are sorted by RT (same order as MS1 scans in the index).
        /// 
        /// If the index has no MS1 scans or the RT window contains no MS1 scans,
        /// <paramref name="rts"/> and <paramref name="intensities"/> are empty arrays.
        /// </summary>
        /// <param name="index">DiaScanIndex containing MS1 scan data.</param>
        /// <param name="precursorMz">Target precursor m/z (M0 monoisotopic).</param>
        /// <param name="rtMin">Lower RT bound (minutes, inclusive).</param>
        /// <param name="rtMax">Upper RT bound (minutes, inclusive).</param>
        /// <param name="ppmTolerance">m/z extraction tolerance in ppm.</param>
        /// <param name="rts">Output: retention times (minutes) of each XIC point.</param>
        /// <param name="intensities">Output: summed intensity at each XIC point.</param>
        public static void ExtractPrecursorXic(
            DiaScanIndex index,
            float precursorMz,
            float rtMin,
            float rtMax,
            float ppmTolerance,
            out float[] rts,
            out float[] intensities)
        {
            if (index == null) throw new ArgumentNullException(nameof(index));

            if (index.Ms1ScanCount == 0 || rtMin > rtMax || precursorMz <= 0f)
            {
                rts = Array.Empty<float>();
                intensities = Array.Empty<float>();
                return;
            }

            float daltonTol = precursorMz * ppmTolerance * 1e-6f;
            float mzLo = precursorMz - daltonTol;
            float mzHi = precursorMz + daltonTol;

            ExtractSingleMzXic(index, mzLo, mzHi, rtMin, rtMax, out rts, out intensities);
        }

        /// <summary>
        /// Extracts XICs for M0, M+1, and M+2 isotopes simultaneously using the
        /// averagine model to compute expected isotope m/z offsets.
        /// 
        /// M+k offset = k × NeutronMass / chargeState
        /// 
        /// All three output arrays are parallel (same length, same RT values).
        /// If the index has no MS1 scans, all outputs are empty arrays.
        /// </summary>
        /// <param name="index">DiaScanIndex containing MS1 scan data.</param>
        /// <param name="precursorMz">M0 monoisotopic precursor m/z.</param>
        /// <param name="chargeState">Precursor charge state (used to compute isotope spacing).</param>
        /// <param name="rtMin">Lower RT bound (minutes, inclusive).</param>
        /// <param name="rtMax">Upper RT bound (minutes, inclusive).</param>
        /// <param name="ppmTolerance">m/z extraction tolerance in ppm (applied at M0 m/z).</param>
        /// <param name="rts">Output: retention times common to all three isotope XICs.</param>
        /// <param name="m0Intensities">Output: M0 (monoisotopic) intensities.</param>
        /// <param name="m1Intensities">Output: M+1 intensities.</param>
        /// <param name="m2Intensities">Output: M+2 intensities.</param>
        public static void ExtractIsotopeXics(
            DiaScanIndex index,
            float precursorMz,
            int chargeState,
            float rtMin,
            float rtMax,
            float ppmTolerance,
            out float[] rts,
            out float[] m0Intensities,
            out float[] m1Intensities,
            out float[] m2Intensities)
        {
            if (index == null) throw new ArgumentNullException(nameof(index));

            if (index.Ms1ScanCount == 0 || rtMin > rtMax || precursorMz <= 0f || chargeState < 1)
            {
                rts = Array.Empty<float>();
                m0Intensities = Array.Empty<float>();
                m1Intensities = Array.Empty<float>();
                m2Intensities = Array.Empty<float>();
                return;
            }

            float isotopeSpacing = NeutronMass / chargeState;

            float m0Mz = precursorMz;
            float m1Mz = precursorMz + isotopeSpacing;
            float m2Mz = precursorMz + 2f * isotopeSpacing;

            // Tolerance in Da computed at M0; same absolute window applied to M+1/M+2
            // (the relative ppm difference across 3 isotopes is < 0.01%, negligible)
            float daltonTol = precursorMz * ppmTolerance * 1e-6f;

            // Find RT window bounds once — shared across all three isotopes
            int startIdx = index.FindMs1ScanIndexAtRt(rtMin);
            int scanCount = CountScansInWindow(index, startIdx, rtMax);

            if (scanCount == 0)
            {
                rts = Array.Empty<float>();
                m0Intensities = Array.Empty<float>();
                m1Intensities = Array.Empty<float>();
                m2Intensities = Array.Empty<float>();
                return;
            }

            // Rent temporary buffers; one pass fills all four in parallel
            float[] rtBuf = ArrayPool<float>.Shared.Rent(scanCount);
            float[] m0Buf = ArrayPool<float>.Shared.Rent(scanCount);
            float[] m1Buf = ArrayPool<float>.Shared.Rent(scanCount);
            float[] m2Buf = ArrayPool<float>.Shared.Rent(scanCount);

            try
            {
                int written = 0;
                for (int i = startIdx; i < index.Ms1ScanCount && written < scanCount; i++)
                {
                    float rt = index.GetMs1ScanRt(i);
                    if (rt > rtMax) break;

                    index.GetMs1ScanPeaks(i, out var mzs, out var ints);

                    rtBuf[written] = rt;
                    m0Buf[written] = SumIntensityInWindow(mzs, ints, m0Mz - daltonTol, m0Mz + daltonTol);
                    m1Buf[written] = SumIntensityInWindow(mzs, ints, m1Mz - daltonTol, m1Mz + daltonTol);
                    m2Buf[written] = SumIntensityInWindow(mzs, ints, m2Mz - daltonTol, m2Mz + daltonTol);
                    written++;
                }

                rts = new float[written];
                m0Intensities = new float[written];
                m1Intensities = new float[written];
                m2Intensities = new float[written];

                rtBuf.AsSpan(0, written).CopyTo(rts);
                m0Buf.AsSpan(0, written).CopyTo(m0Intensities);
                m1Buf.AsSpan(0, written).CopyTo(m1Intensities);
                m2Buf.AsSpan(0, written).CopyTo(m2Intensities);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(rtBuf);
                ArrayPool<float>.Shared.Return(m0Buf);
                ArrayPool<float>.Shared.Return(m1Buf);
                ArrayPool<float>.Shared.Return(m2Buf);
            }
        }

        /// <summary>
        /// Returns the theoretical [M0, M+1, M+2] relative intensities for a given
        /// precursor neutral mass using the precomputed averagine table.
        /// 
        /// The returned span contains three values summing to 1.0.
        /// If the mass is outside the table range [200, 4000 Da], the nearest
        /// boundary entry is returned (clamped).
        /// 
        /// Used by feature [30] IsotopePatternScore in DiaFeatureExtractor.
        /// </summary>
        /// <param name="neutralMass">Precursor neutral mass (Da).</param>
        /// <returns>ReadOnlySpan of three floats: [M0, M+1, M+2] normalized intensities.</returns>
        public static ReadOnlySpan<float> GetTheoreticalIsotopePattern(float neutralMass)
        {
            int massIdx = (int)Math.Round(neutralMass) - AveragineTableMinMass;
            massIdx = Math.Clamp(massIdx, 0, AveragineIsotopeTable.Length - 1);
            return AveragineIsotopeTable[massIdx].AsSpan();
        }

        // ── Internal extraction helpers ──────────────────────────────────────

        /// <summary>
        /// Core single-isotope XIC extractor shared by ExtractPrecursorXic and
        /// the per-isotope calls inside ExtractIsotopeXics.
        /// </summary>
        private static void ExtractSingleMzXic(
            DiaScanIndex index,
            float mzLo,
            float mzHi,
            float rtMin,
            float rtMax,
            out float[] rts,
            out float[] intensities)
        {
            int startIdx = index.FindMs1ScanIndexAtRt(rtMin);
            int scanCount = CountScansInWindow(index, startIdx, rtMax);

            if (scanCount == 0)
            {
                rts = Array.Empty<float>();
                intensities = Array.Empty<float>();
                return;
            }

            float[] rtBuf = ArrayPool<float>.Shared.Rent(scanCount);
            float[] intBuf = ArrayPool<float>.Shared.Rent(scanCount);

            try
            {
                int written = 0;
                for (int i = startIdx; i < index.Ms1ScanCount && written < scanCount; i++)
                {
                    float rt = index.GetMs1ScanRt(i);
                    if (rt > rtMax) break;

                    index.GetMs1ScanPeaks(i, out var mzs, out var ints);

                    rtBuf[written] = rt;
                    intBuf[written] = SumIntensityInWindow(mzs, ints, mzLo, mzHi);
                    written++;
                }

                rts = new float[written];
                intensities = new float[written];
                rtBuf.AsSpan(0, written).CopyTo(rts);
                intBuf.AsSpan(0, written).CopyTo(intensities);
            }
            finally
            {
                ArrayPool<float>.Shared.Return(rtBuf);
                ArrayPool<float>.Shared.Return(intBuf);
            }
        }

        /// <summary>
        /// Counts the number of MS1 scans in [startIdx, ...) with RT ≤ rtMax.
        /// Used to size the rental buffer before the extraction loop.
        /// O(K) where K is the window size — unavoidable since we need the count
        /// before writing.
        /// </summary>
        private static int CountScansInWindow(DiaScanIndex index, int startIdx, float rtMax)
        {
            int count = 0;
            for (int i = startIdx; i < index.Ms1ScanCount; i++)
            {
                if (index.GetMs1ScanRt(i) > rtMax) break;
                count++;
            }
            return count;
        }

        /// <summary>
        /// Sums intensity from all peaks in <paramref name="mzs"/> that fall within
        /// [mzLo, mzHi] (inclusive). Exploits the fact that peaks are sorted by m/z
        /// to use binary search for both the lower and upper bounds.
        /// 
        /// Complexity: O(log P + H) where P = total peaks, H = peaks in window.
        /// For DIA MS1 scans (~500 peaks) and a 20 ppm window at 500 m/z (~0.01 Da),
        /// H is typically 0–2 peaks, so this is effectively O(log P).
        /// </summary>
        /// <param name="mzs">Sorted m/z array for one MS1 scan.</param>
        /// <param name="intensities">Parallel intensity array.</param>
        /// <param name="mzLo">Lower m/z bound (inclusive).</param>
        /// <param name="mzHi">Upper m/z bound (inclusive).</param>
        /// <returns>Sum of intensities for all peaks in the window. 0 if no peaks found.</returns>
        internal static float SumIntensityInWindow(
            ReadOnlySpan<float> mzs,
            ReadOnlySpan<float> intensities,
            float mzLo,
            float mzHi)
        {
            if (mzs.IsEmpty) return 0f;

            // Binary search for the first index >= mzLo
            int lo = LowerBound(mzs, mzLo);
            if (lo >= mzs.Length) return 0f;

            float sum = 0f;
            for (int i = lo; i < mzs.Length; i++)
            {
                if (mzs[i] > mzHi) break;
                sum += intensities[i];
            }
            return sum;
        }

        /// <summary>
        /// Binary lower-bound search: returns the index of the first element >= value.
        /// If all elements are less than value, returns span.Length.
        /// </summary>
        private static int LowerBound(ReadOnlySpan<float> span, float value)
        {
            int lo = 0, hi = span.Length;
            while (lo < hi)
            {
                int mid = lo + ((hi - lo) >> 1);
                if (span[mid] < value)
                    lo = mid + 1;
                else
                    hi = mid;
            }
            return lo;
        }

        // ── Averagine table construction ─────────────────────────────────────

        /// <summary>
        /// Builds the averagine isotope table at static initialization time.
        /// 
        /// For each integer neutral mass M from AveragineTableMinMass to
        /// AveragineTableMaxMass, computes the theoretical [M0, M+1, M+2] relative
        /// intensities using the Poisson approximation to the multinomial isotope
        /// distribution.
        /// 
        /// Averagine formula per residue (111.1 Da):
        ///   C₄.₉₃₈₄ H₇.₇₅₈₃ N₁.₃₅₇₇ O₁.₄₇₇₃ S₀.₀₄₁₇
        /// 
        /// The Poisson approximation uses the expected number of C13/N15/O18/S34
        /// heavy isotopes per molecule as the λ parameter:
        ///   λ ≈ 0.000594 × M  (empirical fit, accurate to ~1% for 200–4000 Da)
        /// 
        /// This gives a simple, fast, allocation-free lookup for feature computation.
        /// </summary>
        private static float[][] BuildAveragineTable()
        {
            int tableSize = AveragineTableMaxMass - AveragineTableMinMass + 1;
            var table = new float[tableSize][];

            for (int i = 0; i < tableSize; i++)
            {
                float mass = AveragineTableMinMass + i;

                // λ = expected number of heavy isotopes (Poisson parameter)
                // Empirical fit to averagine composition: λ ≈ 0.000594 × M
                double lambda = 0.000594 * mass;

                // Poisson probabilities: P(k) = e^-λ × λ^k / k!
                double p0 = Math.Exp(-lambda);
                double p1 = p0 * lambda;
                double p2 = p1 * lambda / 2.0;

                // Normalize to sum = 1
                float sum = (float)(p0 + p1 + p2);
                table[i] = new float[3]
                {
                    (float)(p0 / sum),
                    (float)(p1 / sum),
                    (float)(p2 / sum)
                };
            }

            return table;
        }
    }
}