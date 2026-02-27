// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using MassSpectrometry.Dia;
using NUnit.Framework;
using System;
using System.Numerics;

namespace Test.Dia
{
    /// <summary>
    /// Tests for NormalizedDotProductScorer and SpectralAngleScorer.
    /// 
    /// Both scorers are tested against known mathematical results to verify
    /// correctness. Tests cover: perfect match, orthogonal, scale invariance,
    /// partial overlap, edge cases, and SIMD consistency.
    /// </summary>
    [TestFixture]
    public class ScorerTests
    {
        private NormalizedDotProductScorer _dotProduct;
        private SpectralAngleScorer _spectralAngle;

        [SetUp]
        public void Setup()
        {
            _dotProduct = new NormalizedDotProductScorer();
            _spectralAngle = new SpectralAngleScorer();
        }

        // ════════════════════════════════════════════════════════════════════
        //  Normalized Dot Product Tests
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Identical vectors should produce a perfect score of 1.0.
        /// This is the most basic correctness check — cosine of 0° = 1.
        /// </summary>
        [Test]
        public void DotProduct_IdenticalVectors_ReturnsOne()
        {
            float[] a = { 100f, 200f, 300f, 50f };
            float score = _dotProduct.Score(a, a);
            Assert.That(score, Is.EqualTo(1.0f).Within(1e-5f));
        }

        /// <summary>
        /// Orthogonal vectors (no overlap in nonzero positions) should score 0.0.
        /// This verifies the denominator normalization is correct.
        /// </summary>
        [Test]
        public void DotProduct_OrthogonalVectors_ReturnsZero()
        {
            float[] a = { 100f, 0f, 0f };
            float[] b = { 0f, 200f, 0f };
            float score = _dotProduct.Score(a, b);
            Assert.That(score, Is.EqualTo(0.0f).Within(1e-5f));
        }

        /// <summary>
        /// Cosine similarity is scale-invariant: multiplying all intensities by a
        /// constant should not change the score. This is critical because DIA
        /// observed intensities can vary by orders of magnitude from predictions.
        /// </summary>
        [Test]
        public void DotProduct_ScaledVectors_ReturnsOne()
        {
            float[] a = { 100f, 200f, 300f };
            float[] b = { 1000f, 2000f, 3000f }; // 10× scaled
            float score = _dotProduct.Score(a, b);
            Assert.That(score, Is.EqualTo(1.0f).Within(1e-5f));
        }

        /// <summary>
        /// Two vectors with partial overlap should produce an intermediate score.
        /// Known result: a = [1, 0], b = [1, 1] → cos(45°) = 1/√2 ≈ 0.7071
        /// </summary>
        [Test]
        public void DotProduct_PartialOverlap_ReturnsExpectedValue()
        {
            float[] a = { 1f, 0f };
            float[] b = { 1f, 1f };
            float score = _dotProduct.Score(a, b);
            float expected = 1f / MathF.Sqrt(2f); // cos(45°) ≈ 0.7071
            Assert.That(score, Is.EqualTo(expected).Within(1e-4f));
        }

        /// <summary>
        /// A zero vector (all intensities zero) should return 0, not NaN or infinity.
        /// This can occur when no fragment ions are detected.
        /// </summary>
        [Test]
        public void DotProduct_ZeroVector_ReturnsZero()
        {
            float[] a = { 0f, 0f, 0f };
            float[] b = { 100f, 200f, 300f };
            Assert.That(_dotProduct.Score(a, b), Is.EqualTo(0f));
            Assert.That(_dotProduct.Score(b, a), Is.EqualTo(0f));
        }

        /// <summary>
        /// Empty vectors should return 0.
        /// </summary>
        [Test]
        public void DotProduct_EmptyVectors_ReturnsZero()
        {
            Assert.That(_dotProduct.Score(ReadOnlySpan<float>.Empty, ReadOnlySpan<float>.Empty), Is.EqualTo(0f));
        }

        /// <summary>
        /// Verifies that SIMD and scalar paths produce the same result by testing with
        /// a vector longer than Vector&lt;float&gt;.Count (which triggers SIMD) and comparing
        /// against a manually computed result.
        /// </summary>
        [Test]
        public void DotProduct_LongVector_SimdConsistency()
        {
            // Create vectors longer than SIMD width to ensure the SIMD path is exercised
            int length = Vector<float>.Count * 4 + 3; // Ensures both SIMD and scalar tail
            float[] a = new float[length];
            float[] b = new float[length];

            // Known pattern: a[i] = i+1, b[i] = 1 for all i
            // dot = sum(i+1) = n(n+1)/2, normA = sqrt(sum((i+1)^2)), normB = sqrt(n)
            for (int i = 0; i < length; i++)
            {
                a[i] = i + 1;
                b[i] = 1f;
            }

            float score = _dotProduct.Score(a, b);

            // Compute expected analytically
            float dot = length * (length + 1f) / 2f;
            float normA = 0f;
            for (int i = 0; i < length; i++) normA += (i + 1f) * (i + 1f);
            normA = MathF.Sqrt(normA);
            float normB = MathF.Sqrt(length);
            float expected = dot / (normA * normB);

            Assert.That(score, Is.EqualTo(expected).Within(1e-4f));
        }

        // ════════════════════════════════════════════════════════════════════
        //  Spectral Angle Tests
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Identical vectors → 0° angle → SA = 1.0
        /// </summary>
        [Test]
        public void SpectralAngle_IdenticalVectors_ReturnsOne()
        {
            float[] a = { 100f, 200f, 300f, 50f };
            float score = _spectralAngle.Score(a, a);
            Assert.That(score, Is.EqualTo(1.0f).Within(1e-5f));
        }

        /// <summary>
        /// Orthogonal vectors → 90° angle → SA = 0.0
        /// </summary>
        [Test]
        public void SpectralAngle_OrthogonalVectors_ReturnsZero()
        {
            float[] a = { 100f, 0f, 0f };
            float[] b = { 0f, 200f, 0f };
            float score = _spectralAngle.Score(a, b);
            Assert.That(score, Is.EqualTo(0.0f).Within(1e-5f));
        }

        /// <summary>
        /// Scale invariance: same shape at different scale → SA = 1.0
        /// </summary>
        [Test]
        public void SpectralAngle_ScaledVectors_ReturnsOne()
        {
            float[] a = { 100f, 200f, 300f };
            float[] b = { 500f, 1000f, 1500f };
            float score = _spectralAngle.Score(a, b);
            Assert.That(score, Is.EqualTo(1.0f).Within(1e-5f));
        }

        /// <summary>
        /// Partial overlap: a = [1, 0], b = [1, 1] → cos = 1/√2 → angle = 45°
        /// SA = 1 - (2×45°/180°) = 1 - 0.5 = 0.5
        /// </summary>
        [Test]
        public void SpectralAngle_PartialOverlap_ReturnsExpectedValue()
        {
            float[] a = { 1f, 0f };
            float[] b = { 1f, 1f };
            float score = _spectralAngle.Score(a, b);
            Assert.That(score, Is.EqualTo(0.5f).Within(1e-4f));
        }

        /// <summary>
        /// Zero and empty vectors return 0.
        /// </summary>
        [Test]
        public void SpectralAngle_ZeroAndEmpty_ReturnsZero()
        {
            float[] a = { 0f, 0f };
            float[] b = { 100f, 200f };
            Assert.That(_spectralAngle.Score(a, b), Is.EqualTo(0f));
            Assert.That(_spectralAngle.Score(ReadOnlySpan<float>.Empty, ReadOnlySpan<float>.Empty), Is.EqualTo(0f));
        }

        // ════════════════════════════════════════════════════════════════════
        //  Cross-scorer consistency
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Both scorers should agree on the relative ordering of spectrum pairs.
        /// If pair A is more similar than pair B under dot product, it should also
        /// be more similar under spectral angle. This is a key property for any
        /// downstream ranking or FDR estimation.
        /// </summary>
        [Test]
        public void BothScorers_PreserveRelativeOrdering()
        {
            // Pair A: very similar
            float[] obs = { 100f, 200f, 300f, 400f };
            float[] goodMatch = { 110f, 190f, 310f, 390f };

            // Pair B: less similar
            float[] poorMatch = { 400f, 100f, 200f, 300f }; // shuffled intensities

            float dpGood = _dotProduct.Score(obs, goodMatch);
            float dpPoor = _dotProduct.Score(obs, poorMatch);
            float saGood = _spectralAngle.Score(obs, goodMatch);
            float saPoor = _spectralAngle.Score(obs, poorMatch);

            Assert.That(dpGood, Is.GreaterThan(dpPoor), "Dot product should rank good match higher");
            Assert.That(saGood, Is.GreaterThan(saPoor), "Spectral angle should rank good match higher");
        }
    }
}
