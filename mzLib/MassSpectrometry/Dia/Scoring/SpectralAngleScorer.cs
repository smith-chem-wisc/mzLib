// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Computes the spectral contrast angle (SA) between observed and expected
    /// fragment intensity vectors.
    /// 
    /// The spectral angle is defined as:
    ///   SA = 1 - (2 × arccos(cosine_similarity) / π)
    /// 
    /// This transforms the cosine similarity into an angular measure where:
    ///   1.0 = identical (0° angle between vectors)
    ///   0.0 = orthogonal (90° angle)
    /// 
    /// Compared to the normalized dot product:
    ///   - More sensitive to small differences between similar spectra
    ///   - Better discrimination in the high-similarity range (0.9–1.0)
    ///   - Used by DIA-NN and other modern DIA tools
    /// 
    /// This scorer delegates the cosine computation to the same SIMD-accelerated
    /// logic as NormalizedDotProductScorer, then applies the angular transform.
    /// </summary>
    public sealed class SpectralAngleScorer : IScorer
    {
        /// <summary>
        /// Computes the spectral contrast angle between observed and expected vectors.
        /// Returns a value in [0, 1] where 1 is perfect match.
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

                dotProduct = VectorSum(vDot);
                normObs = VectorSum(vNormObs);
                normExp = VectorSum(vNormExp);
            }

            for (int i = simdEnd; i < length; i++)
            {
                float o = observed[i];
                float e = expected[i];
                dotProduct += o * e;
                normObs += o * o;
                normExp += e * e;
            }

            if (normObs <= 0f || normExp <= 0f) return 0f;

            float cosineSim = dotProduct / (MathF.Sqrt(normObs) * MathF.Sqrt(normExp));
            cosineSim = Math.Clamp(cosineSim, -1f, 1f); // Clamp for arccos safety

            // SA = 1 - (2 * arccos(cos) / π)
            float angle = MathF.Acos(cosineSim);
            float score = 1f - (2f * angle / MathF.PI);

            return Math.Clamp(score, 0f, 1f);
        }

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
