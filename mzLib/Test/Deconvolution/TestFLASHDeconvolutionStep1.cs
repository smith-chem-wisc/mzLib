using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Critical unit tests for FLASHDeconv Step 1:
    ///   - GetLogMz (log transformation)
    ///   - BuildLogMzPeaks (spectrum transformation)
    ///   - BuildUniversalPattern (charge pattern construction)
    ///   - BuildHarmonicPatterns (harmonic artifact patterns)
    ///   - SpectralSimilarity.CosineOfAlignedVectors (new static method)
    ///
    /// Mathematical invariants verified here are the foundation on which
    /// all subsequent steps depend. If these fail, nothing above them is valid.
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvolutionStep1
    {
        // ══════════════════════════════════════════════════════════════════════
        // 1. GetLogMz — the fundamental transformation
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void GetLogMz_PositiveMode_ReturnsLogOfUnchargedMz()
        {
            // For a peak at m/z = mz in positive mode:
            // uncharged = mz - proton
            // logMz     = log(mz - proton)
            double mz = 500.0;
            double expected = Math.Log(mz - Constants.ProtonMass);

            double result = FLASHDeconvolutionAlgorithm.GetLogMz(mz, Polarity.Positive);

            Assert.That(result, Is.EqualTo(expected).Within(1e-12));
        }

        [Test]
        public void GetLogMz_NegativeMode_AddsProtonBeforeLog()
        {
            // In negative mode the charge carrier is removed (proton subtracted),
            // so uncharged = mz + proton → logMz = log(mz + proton)
            double mz = 500.0;
            double expected = Math.Log(mz + Constants.ProtonMass);

            double result = FLASHDeconvolutionAlgorithm.GetLogMz(mz, Polarity.Negative);

            Assert.That(result, Is.EqualTo(expected).Within(1e-12));
        }

        [Test]
        public void GetLogMz_UnphysicalMz_ReturnsNegativeInfinity()
        {
            // An m/z smaller than the proton mass is unphysical in positive mode
            double unphysicalMz = Constants.ProtonMass * 0.5;

            double result = FLASHDeconvolutionAlgorithm.GetLogMz(unphysicalMz, Polarity.Positive);

            Assert.That(result, Is.EqualTo(double.NegativeInfinity));
        }

        /// <summary>
        /// Key invariant: two peaks from the SAME protein at charge states c1 and c2
        /// must differ in log-m/z by exactly log(c2) - log(c1), independent of mass.
        /// This is the whole basis of the universal pattern.
        /// </summary>
        [Test]
        [TestCase(10000.0, 5, 10)]   // 10 kDa protein, charges 5 and 10
        [TestCase(50000.0, 10, 20)]  // 50 kDa protein, charges 10 and 20
        [TestCase(12223.0, 8, 9)]    // Cytochrome C mass range, charges 8 and 9
        public void GetLogMz_TwoChargesOfSameMass_DifferByLogRatio(
            double mass, int z1, int z2)
        {
            double mz1 = mass.ToMz(z1);
            double mz2 = mass.ToMz(z2);

            double logMz1 = FLASHDeconvolutionAlgorithm.GetLogMz(mz1, Polarity.Positive);
            double logMz2 = FLASHDeconvolutionAlgorithm.GetLogMz(mz2, Polarity.Positive);

            double expected = Math.Log(z2) - Math.Log(z1);
            double actual = logMz1 - logMz2;  // higher charge → smaller logMz

            Assert.That(actual, Is.EqualTo(expected).Within(1e-10),
                $"Mass={mass}, z1={z1}, z2={z2}: log-m/z difference should equal log(z2/z1)");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. BuildLogMzPeaks — spectrum transformation
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void BuildLogMzPeaks_BasicSpectrum_ReturnsSortedByLogMz()
        {
            // Arrange: deliberately give mz values out of ascending logMz order
            // (they ARE ascending in mz, but we verify the sort is on logMz not mz)
            var spectrum = new MzSpectrum(
                new[] { 300.0, 500.0, 800.0 },
                new[] { 1000.0, 2000.0, 500.0 },
                shouldCopy: false);
            var range = spectrum.Range;

            var peaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(spectrum, range, Polarity.Positive);

            Assert.That(peaks.Count, Is.EqualTo(3));
            for (int i = 1; i < peaks.Count; i++)
                Assert.That(peaks[i].LogMz, Is.GreaterThanOrEqualTo(peaks[i - 1].LogMz));
        }

        [Test]
        public void BuildLogMzPeaks_ZeroIntensityPeaks_AreExcluded()
        {
            var spectrum = new MzSpectrum(
                new[] { 300.0, 500.0, 800.0 },
                new[] { 0.0, 2000.0, 0.0 },
                shouldCopy: false);

            var peaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                spectrum, spectrum.Range, Polarity.Positive);

            Assert.That(peaks.Count, Is.EqualTo(1));
            Assert.That(peaks[0].Mz, Is.EqualTo(500.0).Within(1e-10));
        }

        [Test]
        public void BuildLogMzPeaks_PeaksOutsideRange_AreExcluded()
        {
            var spectrum = new MzSpectrum(
                new[] { 300.0, 500.0, 800.0 },
                new[] { 1000.0, 2000.0, 500.0 },
                shouldCopy: false);
            var range = new MzRange(400.0, 600.0); // only 500 is inside

            var peaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                spectrum, range, Polarity.Positive);

            Assert.That(peaks.Count, Is.EqualTo(1));
            Assert.That(peaks[0].Mz, Is.EqualTo(500.0).Within(1e-10));
        }

        [Test]
        public void BuildLogMzPeaks_EmptySpectrum_ReturnsEmptyList()
        {
            var spectrum = new MzSpectrum(
                Array.Empty<double>(), Array.Empty<double>(), shouldCopy: false);

            var peaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                spectrum, new MzRange(0, 2000), Polarity.Positive);

            Assert.That(peaks, Is.Empty);
        }

        [Test]
        public void BuildLogMzPeaks_LogMzValues_MatchGetLogMz()
        {
            // Every LogMzPeak.LogMz must equal GetLogMz(mz, polarity) exactly
            var mzs = new[] { 400.0, 600.0, 900.0 };
            var intensities = new[] { 500.0, 1000.0, 250.0 };
            var spectrum = new MzSpectrum(mzs, intensities, shouldCopy: false);

            var peaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                spectrum, spectrum.Range, Polarity.Positive);

            foreach (var peak in peaks)
            {
                double expected = FLASHDeconvolutionAlgorithm.GetLogMz(peak.Mz, Polarity.Positive);
                Assert.That(peak.LogMz, Is.EqualTo(expected).Within(1e-12));
            }
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. BuildUniversalPattern — mathematical invariants
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void BuildUniversalPattern_Length_EqualsChargeRange()
        {
            var p = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbsCharge: 1, chargeRange: 60);
            Assert.That(p.Length, Is.EqualTo(60));
        }

        [Test]
        public void BuildUniversalPattern_Values_AreNegativeLogOfCharge()
        {
            // U[j] = −log(minCharge + j)
            int minCharge = 2;
            int chargeRange = 10;
            var pattern = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minCharge, chargeRange);

            for (int j = 0; j < chargeRange; j++)
            {
                double expected = -Math.Log(minCharge + j);
                Assert.That(pattern[j], Is.EqualTo(expected).Within(1e-12),
                    $"pattern[{j}] should be −log({minCharge + j})");
            }
        }

        [Test]
        public void BuildUniversalPattern_IsStrictlyDecreasing()
        {
            // −log(c) decreases as c increases, so the pattern should be strictly decreasing
            var pattern = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(1, 50);
            for (int i = 1; i < pattern.Length; i++)
                Assert.That(pattern[i], Is.LessThan(pattern[i - 1]),
                    $"pattern[{i}] should be less than pattern[{i - 1}]");
        }

        /// <summary>
        /// The universal pattern spacing equals the log ratio of consecutive charges.
        /// This is the mathematical guarantee that makes FLASHDeconv work.
        /// </summary>
        [Test]
        public void BuildUniversalPattern_ConsecutiveSpacing_EqualsLogRatioOfCharges()
        {
            int minCharge = 1;
            var pattern = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minCharge, 20);

            for (int j = 1; j < pattern.Length; j++)
            {
                int c1 = minCharge + j - 1;
                int c2 = minCharge + j;
                double expectedSpacing = Math.Log(c1) - Math.Log(c2); // negative
                double actualSpacing = pattern[j] - pattern[j - 1];

                Assert.That(actualSpacing, Is.EqualTo(expectedSpacing).Within(1e-12),
                    $"Spacing at j={j} (charges {c1}→{c2}) should equal log({c1}/{c2})");
            }
        }

        [Test]
        public void BuildUniversalPattern_InvalidChargeRange_Throws()
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => FLASHDeconvolutionAlgorithm.BuildUniversalPattern(1, 0));
        }

        [Test]
        public void BuildUniversalPattern_InvalidMinCharge_Throws()
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => FLASHDeconvolutionAlgorithm.BuildUniversalPattern(0, 10));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 4. BuildHarmonicPatterns — shape and relationship to universal pattern
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void BuildHarmonicPatterns_Shape_IsCorrect()
        {
            int chargeRange = 30;
            var hp = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(1, chargeRange);

            // Should have one row per harmonic denominator {2, 3, 5}
            Assert.That(hp.Length, Is.EqualTo(3));
            foreach (var row in hp)
                Assert.That(row.Length, Is.EqualTo(chargeRange));
        }

        [Test]
        public void BuildHarmonicPatterns_Values_LieBetweenConsecutiveUniversalPatternEntries()
        {
            // Each harmonic pattern value for charge index j should fall BETWEEN
            // universalPattern[j-1] and universalPattern[j] (in log-m/z space).
            // This is what "interleaved" means — the harmonic artifact peaks sit
            // between consecutive true-charge peaks.
            int minCharge = 1;
            int chargeRange = 20;
            var up = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minCharge, chargeRange);
            var hp = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minCharge, chargeRange);

            // Start from j=1 because j=0 has no previous entry
            for (int k = 0; k < hp.Length; k++)
            {
                for (int j = 1; j < chargeRange; j++)
                {
                    double lower = up[j];     // more negative (higher charge)
                    double upper = up[j - 1]; // less negative (lower charge)

                    // −log values are all negative; up[j] < up[j-1]
                    // harmonic value should be strictly between them
                    Assert.That(hp[k][j], Is.GreaterThan(lower).And.LessThan(upper),
                        $"Harmonic[{k}][{j}] = {hp[k][j]:F6} should be in " +
                        $"({lower:F6}, {upper:F6})");
                }
            }
        }

        [Test]
        public void BuildHarmonicPatterns_AllValuesAreFinite()
        {
            var hp = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(1, 60);
            foreach (var row in hp)
                foreach (var v in row)
                    Assert.That(double.IsFinite(v), Is.True,
                        $"Harmonic value {v} should be finite");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 5. SpectralSimilarity.CosineOfAlignedVectors — new static method
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void CosineOfAlignedVectors_IdenticalVectors_ReturnsOne()
        {
            double[] a = { 1.0, 2.0, 3.0, 4.0 };
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(a, a),
                Is.EqualTo(1.0).Within(1e-12));
        }

        [Test]
        public void CosineOfAlignedVectors_OrthogonalVectors_ReturnsZero()
        {
            double[] a = { 1.0, 0.0, 0.0 };
            double[] b = { 0.0, 1.0, 0.0 };
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(a, b),
                Is.EqualTo(0.0).Within(1e-12));
        }

        [Test]
        public void CosineOfAlignedVectors_ScaledVector_ReturnsOne()
        {
            // Scaling a vector does not change the cosine similarity
            double[] a = { 1.0, 2.0, 3.0 };
            double[] b = { 5.0, 10.0, 15.0 }; // a * 5
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(a, b),
                Is.EqualTo(1.0).Within(1e-12));
        }

        [Test]
        public void CosineOfAlignedVectors_KnownValue_IsCorrect()
        {
            // [1,1,0] · [1,0,1] = 1; norms both √2 → cosine = 1/2 = 0.5
            double[] a = { 1.0, 1.0, 0.0 };
            double[] b = { 1.0, 0.0, 1.0 };
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(a, b),
                Is.EqualTo(0.5).Within(1e-12));
        }

        [Test]
        public void CosineOfAlignedVectors_DifferentLengths_ReturnsZero()
        {
            double[] a = { 1.0, 2.0 };
            double[] b = { 1.0, 2.0, 3.0 };
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(a, b),
                Is.EqualTo(0.0));
        }

        [Test]
        public void CosineOfAlignedVectors_EmptyVectors_ReturnsZero()
        {
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(
                ReadOnlySpan<double>.Empty, ReadOnlySpan<double>.Empty),
                Is.EqualTo(0.0));
        }

        [Test]
        public void CosineOfAlignedVectors_AllZeroVector_ReturnsZero()
        {
            double[] a = { 0.0, 0.0, 0.0 };
            double[] b = { 1.0, 2.0, 3.0 };
            Assert.That(SpectralSimilarity.CosineOfAlignedVectors(a, b),
                Is.EqualTo(0.0));
        }

        [Test]
        public void CosineOfAlignedVectors_ResultIsInZeroToOneRange()
        {
            // For non-negative intensity arrays (as in isotope distributions),
            // cosine similarity should always be in [0, 1]
            double[] observed = { 0.1, 0.5, 1.0, 0.7, 0.3 };
            double[] theoretical = { 0.05, 0.4, 0.9, 0.6, 0.2 };

            double result = SpectralSimilarity.CosineOfAlignedVectors(observed, theoretical);

            Assert.That(result, Is.InRange(0.0, 1.0));
        }

        [Test]
        public void CosineOfAlignedVectors_PlainArraysPassWithoutCast()
        {
            // Verifies that double[] implicitly converts to ReadOnlySpan<double>
            // — this is the primary calling pattern in FLASHDeconvolutionAlgorithm
            double[] observed = { 1.0, 2.0, 3.0 };
            double[] theoretical = { 1.0, 2.0, 3.0 };

            // If this compiles and runs, implicit conversion works
            Assert.DoesNotThrow(() =>
                _ = SpectralSimilarity.CosineOfAlignedVectors(observed, theoretical));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 6. Integration: log transform preserves the physics of a real protein
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// For Cytochrome C (monoisotopic ~12223.2 Da), peaks at charge states 9–13
        /// should all map to approximately the same log(mass) value after the transform,
        /// and should match the universal pattern offset at log(mass).
        /// This validates that the math actually works for a real-world case.
        /// </summary>
        [Test]
        public void LogTransform_CytoCPeaks_AllMapToSameLogMass()
        {
            const double cytoMass = 12223.2; // approximate monoisotopic mass Da
            int[] charges = { 9, 10, 11, 12, 13 };
            double tolerance = 1e-6; // log units — very tight

            // Compute the expected log(mass) from the first charge
            double mz0 = cytoMass.ToMz(charges[0]);
            double logMz0 = FLASHDeconvolutionAlgorithm.GetLogMz(mz0, Polarity.Positive);
            double expectedLogMass = logMz0 + Math.Log(charges[0]); // = log(mass)

            foreach (int z in charges)
            {
                double mz = cytoMass.ToMz(z);
                double logMz = FLASHDeconvolutionAlgorithm.GetLogMz(mz, Polarity.Positive);
                double logMass = logMz + Math.Log(z); // recover log(mass)

                Assert.That(logMass, Is.EqualTo(expectedLogMass).Within(tolerance),
                    $"Charge {z}: log(mass) = {logMass} should equal {expectedLogMass}");
            }
        }
    }
}