// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Computes the normalized dot product (cosine similarity) between observed and expected
    /// fragment intensity vectors.
    /// 
    /// This is the standard scoring metric for DIA spectral matching, equivalent to:
    ///   score = (A · B) / (||A|| × ||B||)
    /// 
    /// where A and B are the observed and expected intensity vectors.
    /// 
    /// The result is in [0, 1] where:
    ///   1.0 = identical relative intensity patterns
    ///   0.0 = orthogonal (no similarity)
    /// 
    /// Performance notes:
    ///   - No heap allocations
    ///   - Uses SIMD (System.Numerics.Vector) when vectors are large enough
    ///   - Falls back to scalar loop for short vectors
    ///   - Handles edge cases: zero vectors, single elements, negative values
    /// </summary>
    public sealed class NormalizedDotProductScorer : IScorer
    {
        /// <summary>
        /// Computes cosine similarity between observed and expected intensity vectors.
        /// Both spans must have the same length.
        /// 
        /// Returns 0.0 if either vector has zero magnitude (all zeros or all negative
        /// values after any future preprocessing).
        /// </summary>
        public float Score(ReadOnlySpan<float> observed, ReadOnlySpan<float> expected)
        {
            if (observed.Length != expected.Length)
                throw new ArgumentException("Observed and expected spans must have the same length.");

            int length = observed.Length;
            if (length == 0) return 0f;

            float dotProduct = 0f;
            float normObs = 0f;
            float normExp = 0f;

            // SIMD path: process Vector<float>.Count elements at a time
            int simdWidth = Vector<float>.Count;
            int simdEnd = length - (length % simdWidth);

            if (simdEnd >= simdWidth)
            {
                Vector<float> vDot = Vector<float>.Zero;
                Vector<float> vNormObs = Vector<float>.Zero;
                Vector<float> vNormExp = Vector<float>.Zero;

                for (int i = 0; i < simdEnd; i += simdWidth)
                {
                    var vO = new Vector<float>(observed.Slice(i, simdWidth));
                    var vE = new Vector<float>(expected.Slice(i, simdWidth));

                    vDot += vO * vE;
                    vNormObs += vO * vO;
                    vNormExp += vE * vE;
                }

                // Horizontal sum of SIMD accumulators
                dotProduct = VectorSum(vDot);
                normObs = VectorSum(vNormObs);
                normExp = VectorSum(vNormExp);
            }

            // Scalar tail: handle remaining elements
            for (int i = simdEnd; i < length; i++)
            {
                float o = observed[i];
                float e = expected[i];
                dotProduct += o * e;
                normObs += o * o;
                normExp += e * e;
            }

            // Guard against zero-magnitude vectors
            if (normObs <= 0f || normExp <= 0f) return 0f;

            float denominator = MathF.Sqrt(normObs) * MathF.Sqrt(normExp);
            float score = dotProduct / denominator;

            // Clamp to [0, 1] to handle floating-point rounding
            return Math.Clamp(score, 0f, 1f);
        }

        /// <summary>
        /// Computes the horizontal sum of all elements in a SIMD vector.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static float VectorSum(Vector<float> v)
        {
            float sum = 0f;
            for (int i = 0; i < Vector<float>.Count; i++)
            {
                sum += v[i];
            }
            return sum;
        }
    }
}
