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
    /// Critical unit tests for FLASHDeconv Step 2: candidate mass finding.
    ///
    /// Tests are organised into three groups:
    ///   1. Bin arithmetic (BinNumber / BinValue) — the primitive operations
    ///      on which the entire binning system depends.
    ///   2. FindCandidateMasses with synthetic spectra — controlled inputs where
    ///      we know exactly which masses should be found.
    ///   3. Real-world sanity tests — realistic multi-charge spectra derived
    ///      from known protein masses.
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvolutionStep2
    {
        // ── shared parameters ────────────────────────────────────────────────
        private static FLASHDeconvolutionParameters DefaultParams() =>
            new FLASHDeconvolutionParameters(
                minCharge: 1, maxCharge: 60,
                deconvolutionTolerancePpm: 10.0,
                minIsotopicPeakCount: 3,
                minMassRange: 50.0,
                maxMassRange: 100_000.0);

        // ══════════════════════════════════════════════════════════════════════
        // 1. Bin arithmetic
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void BinNumber_MinValue_ReturnsZero()
        {
            double minValue = 5.0;
            double binMulFactor = 100_000.0;
            Assert.That(FLASHDeconvolutionAlgorithm.BinNumber(minValue, minValue, binMulFactor),
                Is.EqualTo(0));
        }

        [Test]
        public void BinNumber_RoundTripsWithBinValue()
        {
            double minValue = Math.Log(50.0);
            double binMulFactor = 100_000.0;
            double value = Math.Log(12000.0);

            int bin = FLASHDeconvolutionAlgorithm.BinNumber(value, minValue, binMulFactor);
            double recovered = FLASHDeconvolutionAlgorithm.BinValue(bin, minValue, binMulFactor);

            // Round-trip tolerance: ≤ 0.5 / binMulFactor in log space
            Assert.That(Math.Abs(recovered - value), Is.LessThan(0.5 / binMulFactor + 1e-12));
        }

        [Test]
        public void BinNumber_BelowMinValue_ReturnsZero()
        {
            // Values below the minimum should clamp to bin 0, not go negative
            Assert.That(FLASHDeconvolutionAlgorithm.BinNumber(4.0, 5.0, 100_000.0),
                Is.EqualTo(0));
        }

        [Test]
        public void BinNumber_LargerValueGivesLargerBin()
        {
            double minValue = Math.Log(100.0);
            double binMulFactor = 100_000.0;

            int bin1 = FLASHDeconvolutionAlgorithm.BinNumber(Math.Log(1_000.0), minValue, binMulFactor);
            int bin2 = FLASHDeconvolutionAlgorithm.BinNumber(Math.Log(10_000.0), minValue, binMulFactor);

            Assert.That(bin2, Is.GreaterThan(bin1));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. FindCandidateMasses — synthetic spectra
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Helper: build a synthetic spectrum for a single protein mass
        /// at a range of charge states, with no noise.
        /// </summary>
        private static MzSpectrum MakeProteinSpectrum(
            double neutralMass,
            int[] charges,
            double baseIntensity = 10_000.0)
        {
            var mzList = new List<double>();
            var itList = new List<double>();

            foreach (int z in charges)
            {
                mzList.Add(neutralMass.ToMz(z));
                itList.Add(baseIntensity);
            }

            // Sort by mz (MzSpectrum requires ascending order)
            var pairs = mzList.Zip(itList).OrderBy(p => p.First).ToArray();
            return new MzSpectrum(
                pairs.Select(p => p.First).ToArray(),
                pairs.Select(p => p.Second).ToArray(),
                shouldCopy: false);
        }

        [Test]
        public void FindCandidateMasses_SingleProtein_FindsAtLeastOneCandidateNearTrueMass()
        {
            // Cytochrome C: ~12,223 Da, charges 8–14
            const double trueMass = 12_223.2;
            int[] charges = { 8, 9, 10, 11, 12, 13, 14 };
            double tolerancePpm = 10.0;

            var spectrum = MakeProteinSpectrum(trueMass, charges);
            var deconParams = DefaultParams();

            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                      spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);
            double binMulFactor = 1.0 / (tolerancePpm * 1e-6);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, binMulFactor, deconParams);

            // At least one candidate should be within 20 ppm of true mass
            double maxDa = trueMass * 20e-6;
            bool found = candidates.Any(c => Math.Abs(c.Mass - trueMass) <= maxDa);

            Assert.That(found, Is.True,
                $"Expected a candidate within 20 ppm of {trueMass} Da. " +
                $"Candidates: [{string.Join(", ", candidates.Select(c => $"{c.Mass:F1}"))}]");
        }

        [Test]
        public void FindCandidateMasses_EmptySpectrum_ReturnsEmpty()
        {
            var spectrum = new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false);
            var deconParams = DefaultParams();
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                  spectrum, new MzRange(0, 2000), Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, 100_000.0, deconParams);

            Assert.That(candidates, Is.Empty);
        }

        [Test]
        public void FindCandidateMasses_SinglePeak_ReturnsEmpty()
        {
            // A single peak cannot satisfy MinSupportPeakCount = 3
            var spectrum = new MzSpectrum(new[] { 500.0 }, new[] { 1000.0 }, false);
            var deconParams = DefaultParams();
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                  spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, 100_000.0, deconParams);

            Assert.That(candidates, Is.Empty,
                "A single peak should not produce any candidate mass.");
        }

        [Test]
        public void FindCandidateMasses_TwoPeaks_ReturnsEmpty()
        {
            // Two peaks — still below MinSupportPeakCount
            double mass = 5_000.0;
            var spectrum = MakeProteinSpectrum(mass, new[] { 5, 6 });
            var deconParams = DefaultParams();
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, 100_000.0, deconParams);

            Assert.That(candidates, Is.Empty,
                "Two peaks are below the minimum support peak count of 3.");
        }

        [Test]
        public void FindCandidateMasses_ThreePeaks_FindsCandidate()
        {
            // Exactly MinSupportPeakCount = 3 consecutive charge states should be enough
            double mass = 5_000.0;
            var spectrum = MakeProteinSpectrum(mass, new[] { 5, 6, 7 });
            var deconParams = DefaultParams();
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                 spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);
            double binMul = 1.0 / (10.0 * 1e-6);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, binMul, deconParams);

            double maxDa = mass * 20e-6;
            bool found = candidates.Any(c => Math.Abs(c.Mass - mass) <= maxDa);

            Assert.That(found, Is.True,
                $"Three consecutive charge-state peaks should produce a candidate near {mass} Da.");
        }

        [Test]
        public void FindCandidateMasses_AllCandidatesMassesInAllowedRange()
        {
            // All returned candidates must respect MinMassRange and MaxMassRange
            double mass = 12_223.0;
            var spectrum = MakeProteinSpectrum(mass, new[] { 8, 9, 10, 11, 12 });
            var deconParams = DefaultParams(); // MinMassRange=50, MaxMassRange=100,000
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                 spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, 100_000.0, deconParams);

            foreach (var c in candidates)
            {
                Assert.That(c.Mass, Is.GreaterThanOrEqualTo(deconParams.MinMassRange),
                    $"Candidate mass {c.Mass} is below MinMassRange {deconParams.MinMassRange}");
                Assert.That(c.Mass, Is.LessThanOrEqualTo(deconParams.MaxMassRange),
                    $"Candidate mass {c.Mass} is above MaxMassRange {deconParams.MaxMassRange}");
            }
        }

        [Test]
        public void FindCandidateMasses_CandidatesAreSortedAscendingByMass()
        {
            // Two proteins at well-separated masses
            double mass1 = 5_000.0;
            double mass2 = 20_000.0;

            var mzList = new List<double>();
            var itList = new List<double>();
            foreach (int z in new[] { 5, 6, 7 })
            {
                mzList.Add(mass1.ToMz(z));
                itList.Add(10_000.0);
            }
            foreach (int z in new[] { 10, 12, 14 })
            {
                mzList.Add(mass2.ToMz(z));
                itList.Add(10_000.0);
            }
            var pairs = mzList.Zip(itList).OrderBy(p => p.First).ToArray();
            var spectrum = new MzSpectrum(
                pairs.Select(p => p.First).ToArray(),
                pairs.Select(p => p.Second).ToArray(),
                shouldCopy: false);

            var deconParams = DefaultParams();
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                  spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, 100_000.0, deconParams);

            for (int i = 1; i < candidates.Count; i++)
                Assert.That(candidates[i].Mass, Is.GreaterThanOrEqualTo(candidates[i - 1].Mass),
                    "Candidates should be sorted ascending by mass.");
        }

        [Test]
        public void FindCandidateMasses_ChargeRangeIsPlausible()
        {
            // The charge range recorded in the candidate should be a subset of
            // the charges used to generate the spectrum.
            double mass = 12_223.0;
            int[] charges = { 9, 10, 11, 12 };
            var spectrum = MakeProteinSpectrum(mass, charges);
            var deconParams = DefaultParams();
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                                  spectrum, spectrum.Range, Polarity.Positive);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);

            var candidates = FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, 100_000.0, deconParams);

            double maxDa = mass * 30e-6;
            var matching = candidates.Where(c => Math.Abs(c.Mass - mass) <= maxDa).ToList();

            // At least one matching candidate should have charge range within [9,12] ± 2
            bool chargesOk = matching.Any(c =>
                c.MinAbsCharge >= charges.Min() - 2 &&
                c.MaxAbsCharge <= charges.Max() + 2 &&
                c.MinAbsCharge <= c.MaxAbsCharge);

            Assert.That(chargesOk, Is.True,
                $"Candidate charge range should overlap [{charges.Min()},{charges.Max()}]. " +
                $"Got: [{string.Join(", ", matching.Select(c => $"[{c.MinAbsCharge}–{c.MaxAbsCharge}]"))}]");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. Real-world sanity — Deconvoluter integration
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Calls through the full Deconvoluter dispatch path (Step 2 returns empty
        /// because Steps 3–5 are not yet implemented, but it must not throw).
        /// This is the integration guard for the completed infrastructure.
        /// </summary>
        [Test]
        public void Deconvoluter_WithFLASHParams_Step2_DoesNotThrow()
        {
            double mass = 12_223.0;
            var spectrum = MakeProteinSpectrum(mass, new[] { 9, 10, 11, 12, 13 });
            var deconParams = DefaultParams();

            IEnumerable<IsotopicEnvelope> result = null;
            Assert.DoesNotThrow(() =>
                result = Deconvoluter.Deconvolute(spectrum, deconParams).ToList());

            Assert.That(result, Is.Not.Null);
            // Steps 3–5 not yet implemented → still empty; update when implemented
            Assert.That(result, Is.Empty,
                "Steps 3–5 not yet implemented. Update this assertion when they are.");
        }

        [Test]
        public void Deconvoluter_WithFLASHParams_NegativePolarity_DoesNotThrow()
        {
            var spectrum = MakeProteinSpectrum(5_000.0, new[] { 5, 6, 7, 8 });
            var deconParams = new FLASHDeconvolutionParameters(polarity: Polarity.Negative);

            Assert.DoesNotThrow(() =>
                _ = Deconvoluter.Deconvolute(spectrum, deconParams).ToList());
        }
    }
}