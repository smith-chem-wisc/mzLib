using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Tests for FLASHDeconv Steps 3–5: isotope recruitment, cosine scoring, output.
    /// These are the first tests that assert real <see cref="IsotopicEnvelope"/> output
    /// from the algorithm. Step 2 integration tests that asserted empty results must be
    /// updated to reflect that the algorithm now produces output.
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvolutionStep3to5
    {
        // ── helpers ──────────────────────────────────────────────────────────

        private static FLASHDeconvolutionParameters DefaultParams() =>
            new FLASHDeconvolutionParameters(
                minCharge: 1, maxCharge: 60,
                deconvolutionTolerancePpm: 10.0,
                minIsotopicPeakCount: 3,
                minCosineScore: 0.4,  // paper's spectral-level threshold
                minMassRange: 50.0,
                maxMassRange: 100_000.0);

        /// <summary>
        /// Builds a synthetic spectrum with the correct monoisotopic + isotope peaks
        /// for a protein at a given mass and set of charges, using the Averagine model.
        /// Returns both the spectrum and the actual monoisotopic mass used (which is the
        /// lightest peak in the Averagine fine-resolution distribution for the nearest
        /// formula — this is what the algorithm will recover, not the input mass).
        /// </summary>
        private static (MzSpectrum spectrum, double actualMonoMass) MakeSyntheticIsotopeSpectrum(
            double targetMass, int[] charges,
            FLASHDeconvolutionParameters p,
            double baseIntensity = 100_000.0)
        {
            var mzList = new List<double>();
            var itList = new List<double>();

            int avgIdx = p.AverageResidueModel.GetMostIntenseMassIndex(targetMass);
            double[] avgMassesRaw = p.AverageResidueModel.GetAllTheoreticalMasses(avgIdx);
            double[] avgIntensRaw = p.AverageResidueModel.GetAllTheoreticalIntensities(avgIdx);

            // Re-sort to mass-ascending (isotope) order
            var sorted = avgMassesRaw.Zip(avgIntensRaw)
                .OrderBy(pair => pair.First)
                .ToArray();
            double[] isoMasses = sorted.Select(pair => pair.First).ToArray();
            double[] isoIntens = sorted.Select(pair => pair.Second).ToArray();

            // The actual monoisotopic mass is the lightest peak in the Averagine distribution.
            // This equals MostIntenseMasses[avgIdx] - DiffToMonoisotopic[avgIdx].
            // It is NOT necessarily equal to targetMass — Averagine is quantized by formula.
            double actualMonoMass = isoMasses[0];

            foreach (int z in charges)
            {
                for (int n = 0; n < isoMasses.Length; n++)
                {
                    double mz = (actualMonoMass + n * Constants.C13MinusC12).ToMz(z);
                    double intensity = baseIntensity * isoIntens[n];
                    if (intensity < baseIntensity * 0.001) continue;
                    mzList.Add(mz);
                    itList.Add(intensity);
                }
            }

            var pairs = mzList.Zip(itList).OrderBy(pair => pair.First).ToArray();
            var spectrum = new MzSpectrum(
                pairs.Select(pair => pair.First).ToArray(),
                pairs.Select(pair => pair.Second).ToArray(),
                shouldCopy: false);
            return (spectrum, actualMonoMass);
        }

        /// <summary>Convenience overload that discards actualMonoMass — for tests that only need the spectrum.</summary>
        private static MzSpectrum Spectrum(double targetMass, int[] charges, FLASHDeconvolutionParameters p,
            double baseIntensity = 100_000.0)
            => MakeSyntheticIsotopeSpectrum(targetMass, charges, p, baseIntensity).spectrum;

        // ══════════════════════════════════════════════════════════════════════
        // 1. Basic output: algorithm now returns non-empty results
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_SyntheticCytoC_ReturnsAtLeastOneEnvelope()
        {
            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);
            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();
            Assert.That(results, Is.Not.Empty,
                "A perfect synthetic spectrum should produce at least one envelope");
        }

        [Test]
        public void Deconvolute_SyntheticCytoC_FindsMassWithin20Ppm()
        {
            // The algorithm recovers formula.MonoisotopicMass (= apex - DiffToMonoisotopic),
            // which is 12229.16 for the Averagine nearest to 12223.2 — NOT 12223.2 itself.
            // Diagnostic tests A1-C2 confirmed this. We compare against actualMonoMass.
            var deconParams = DefaultParams();
            var (spectrum, actualMono) = MakeSyntheticIsotopeSpectrum(
                12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            double toleranceDa = actualMono * 20e-6;
            bool found = results.Any(e => Math.Abs(e.MonoisotopicMass - actualMono) <= toleranceDa);
            Assert.That(found, Is.True,
                $"Expected an envelope within 20 ppm of actualMono={actualMono:F4} Da. " +
                $"Got: [{string.Join(", ", results.Select(e => $"{e.MonoisotopicMass:F2}"))}]");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. Score is in valid range
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_AllEnvelopes_HaveScoreAboveMinCosine()
        {

            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            foreach (var env in results)
                Assert.That(env.Score, Is.GreaterThanOrEqualTo(deconParams.MinCosineScore),
                    $"Envelope score {env.Score:F4} is below MinCosineScore {deconParams.MinCosineScore}");
        }

        [Test]
        public void Deconvolute_AllEnvelopes_HaveScoreAtMostOne()
        {

            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            foreach (var env in results)
                Assert.That(env.Score, Is.LessThanOrEqualTo(1.0 + 1e-9),
                    $"Cosine score {env.Score:F6} exceeds 1.0");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. MinIsotopicPeakCount filter is respected
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_AllEnvelopes_HaveAtLeastMinIsotopicPeakCount()
        {

            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            foreach (var env in results)
                Assert.That(env.Peaks.Count, Is.GreaterThanOrEqualTo(deconParams.MinIsotopicPeakCount),
                    $"Envelope has {env.Peaks.Count} peaks, below MinIsotopicPeakCount {deconParams.MinIsotopicPeakCount}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 4. Mass range filter is respected
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_AllEnvelopes_MassWithinDeclaredRange()
        {

            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            foreach (var env in results)
            {
                Assert.That(env.MonoisotopicMass, Is.GreaterThanOrEqualTo(deconParams.MinMassRange),
                    $"MonoisotopicMass {env.MonoisotopicMass:F1} is below MinMassRange");
                Assert.That(env.MonoisotopicMass, Is.LessThanOrEqualTo(deconParams.MaxMassRange),
                    $"MonoisotopicMass {env.MonoisotopicMass:F1} is above MaxMassRange");
            }
        }

        // ══════════════════════════════════════════════════════════════════════
        // 5. Charge sign matches polarity
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_PositivePolarity_AllChargesPositive()
        {

            var deconParams = DefaultParams();  // default = Polarity.Positive
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            foreach (var env in results)
                Assert.That(env.Charge, Is.GreaterThan(0),
                    $"Expected positive charge in positive mode, got {env.Charge}");
        }

        [Test]
        public void Deconvolute_NegativePolarity_AllChargesNegative()
        {

            var deconParams = new FLASHDeconvolutionParameters(
                polarity: Polarity.Negative,
                minCosineScore: 0.4);
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            // If any results exist (negative mode spectrum from positive-mode peaks may be empty),
            // all charges must be negative.
            foreach (var env in results)
                Assert.That(env.Charge, Is.LessThan(0),
                    $"Expected negative charge in negative mode, got {env.Charge}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 6. MinCosineScore filter works
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_HighMinCosineScore_ReducesOrEliminatesOutput()
        {

            var paramsLow = new FLASHDeconvolutionParameters(minCosineScore: 0.1);
            var paramsHigh = new FLASHDeconvolutionParameters(minCosineScore: 0.99);

            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, paramsLow);

            var resultsLow = Deconvoluter.Deconvolute(spectrum, paramsLow).ToList();
            var resultsHigh = Deconvoluter.Deconvolute(spectrum, paramsHigh).ToList();

            Assert.That(resultsHigh.Count, Is.LessThanOrEqualTo(resultsLow.Count),
                "A stricter cosine threshold should produce fewer or equal envelopes");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 7. Empty spectrum still returns empty (regression guard)
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_EmptySpectrum_StillReturnsEmpty()
        {
            var spectrum = new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false);
            var deconParams = DefaultParams();

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            Assert.That(results, Is.Empty);
        }

        // ══════════════════════════════════════════════════════════════════════
        // 8. TotalIntensity is positive
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvolute_AllEnvelopes_TotalIntensityIsPositive()
        {

            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            var results = Deconvoluter.Deconvolute(spectrum, deconParams).ToList();

            foreach (var env in results)
                Assert.That(env.TotalIntensity, Is.GreaterThan(0.0),
                    "TotalIntensity must be positive");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 9. Update Step 2 integration test — no longer empty
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvoluter_WithFLASHParams_Step3to5_DoesNotThrowAndReturnsResults()
        {
            // This replaces the Step 2 test that asserted Is.Empty.
            // With Steps 3–5 complete, a realistic spectrum should now yield results.

            var deconParams = DefaultParams();
            var spectrum = Spectrum(12_223.2, new[] { 9, 10, 11, 12, 13 }, deconParams);

            IEnumerable<IsotopicEnvelope> result = null;
            Assert.DoesNotThrow(() =>
                result = Deconvoluter.Deconvolute(spectrum, deconParams).ToList());

            Assert.That(result, Is.Not.Null);
            // With a synthetic perfect-isotope spectrum, we expect results.
            // If this fails, check Steps 3–5 for regressions.
            Assert.That(result.Any(), Is.True,
                "Steps 3–5 implemented: a synthetic isotope spectrum should yield at least one envelope");
        }
    }
}