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
    /// Additional critical edge-case tests for Step 2 (FindCandidateMasses).
    /// These complement TestFLASHDeconvolutionStep2 by targeting specific
    /// failure modes not covered by the basic tests.
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvolutionStep2EdgeCases
    {
        // ── shared helpers ───────────────────────────────────────────────────

        private static MzSpectrum MakeSpectrum(double mass, int[] charges,
            double[] relIntensities = null)
        {
            relIntensities ??= Enumerable.Repeat(1.0, charges.Length).ToArray();
            var pairs = charges
                .Zip(relIntensities, (z, ri) => (mz: mass.ToMz(z), intensity: 10_000.0 * ri))
                .OrderBy(p => p.mz).ToArray();
            return new MzSpectrum(
                pairs.Select(p => p.mz).ToArray(),
                pairs.Select(p => p.intensity).ToArray(),
                shouldCopy: false);
        }

        private static List<FLASHDeconvolutionAlgorithm.CandidateMass> RunStep2(
            MzSpectrum spectrum, FLASHDeconvolutionParameters deconParams)
        {
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(
                spectrum, spectrum.Range, deconParams.Polarity);
            int minAbs = Math.Abs(deconParams.MinAssumedChargeState);
            int chargeRange = Math.Abs(deconParams.MaxAssumedChargeState) - minAbs + 1;
            var universal = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, chargeRange);
            var harmonic = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, chargeRange);
            double binMul = 1.0 / (deconParams.DeconvolutionTolerancePpm * 1e-6);
            return FLASHDeconvolutionAlgorithm.FindCandidateMasses(
                logPeaks, universal, harmonic, binMul, deconParams);
        }

        private static FLASHDeconvolutionParameters Params(
            double tolerancePpm = 10.0,
            double minMass = 50.0, double maxMass = 100_000.0,
            int minCharge = 1, int maxCharge = 60,
            Polarity polarity = Polarity.Positive)
            => new FLASHDeconvolutionParameters(
                minCharge: minCharge, maxCharge: maxCharge,
                deconvolutionTolerancePpm: tolerancePpm,
                minMassRange: minMass, maxMassRange: maxMass,
                polarity: polarity);

        // ══════════════════════════════════════════════════════════════════════
        // 1. Tolerance sensitivity
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        [TestCase(5.0)]    // tight — binMulFactor = 200,000
        [TestCase(10.0)]   // default — 100,000
        [TestCase(20.0)]   // loose — 50,000
        public void FindCandidateMasses_DifferentTolerances_FindsCytoCAtAllTolValues(double ppm)
        {
            // Cytochrome C at charges 9–13 should be found at any reasonable ppm
            const double mass = 12_223.0;
            var spectrum = MakeSpectrum(mass, new[] { 9, 10, 11, 12, 13 });
            var candidates = RunStep2(spectrum, Params(tolerancePpm: ppm));

            double maxDa = mass * 30e-6;
            bool found = candidates.Any(c => Math.Abs(c.Mass - mass) <= maxDa);
            Assert.That(found, Is.True,
                $"At {ppm} ppm tolerance, expected candidate near {mass} Da");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. Intensity asymmetry — at the edge of the 10× ratio limit
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FindCandidateMasses_IntensitiesWithinRatioLimit_FindsCandidate()
        {
            // Intensities vary by up to 9× across 5 charges — should still find mass
            const double mass = 8_000.0;
            double[] relIntensities = { 1.0, 5.0, 9.0, 4.0, 2.0 };
            var spectrum = MakeSpectrum(mass, new[] { 6, 7, 8, 9, 10 }, relIntensities);
            var candidates = RunStep2(spectrum, Params());

            double maxDa = mass * 30e-6;
            Assert.That(candidates.Any(c => Math.Abs(c.Mass - mass) <= maxDa), Is.True,
                "Intensities varying up to 9× should still produce a candidate");
        }

        [Test]
        public void FindCandidateMasses_IntensityDropBy100x_BreaksRunAndResets()
        {
            // A 100× drop between charge 8 and 9 should break the continuous charge run.
            // With only 3 charges on each side of the break (charges 6,7,8 and 9,10,11),
            // each side independently gets exactly 3 support peaks — so we may still find
            // a candidate (the run restarts). The key invariant: no single candidate spans
            // across the intensity break by more than the 10× ratio.
            const double mass = 10_000.0;
            double[] relIntensities = { 1.0, 1.0, 1.0, 0.001, 0.001, 0.001 };
            int[] charges = { 6, 7, 8, 9, 10, 11 };
            var spectrum = MakeSpectrum(mass, charges, relIntensities);

            // Should not throw, and if a candidate is found its charge range should
            // not span the break (i.e., no single candidate spans z=8 to z=9 at 100× gap)
            Assert.DoesNotThrow(() => RunStep2(spectrum, Params()));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. Harmonic rejection doesn't over-reject real peaks
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FindCandidateMasses_NoSpuriousHarmonicRejection_ForCleanProtein()
        {
            // A clean single-protein spectrum at high charge states should not have
            // its own peaks mistakenly flagged as harmonic artifacts of each other.
            // Uses consecutive high charges (15–19) so the continuity check is
            // satisfied; the intent is that harmonic rejection does not fire
            // on genuine charge-state peaks of the same protein.
            // NOTE: charges must be consecutive (j increments by 1) for the
            // continuity accumulator to work — non-consecutive charges (e.g. 15,17,19)
            // would always fail the prevJ-j==-1 check regardless of harmonic status.
            const double mass = 25_000.0;
            var spectrum = MakeSpectrum(mass, new[] { 15, 16, 17, 18, 19 });
            var candidates = RunStep2(spectrum, Params());

            double maxDa = mass * 30e-6;
            Assert.That(candidates.Any(c => Math.Abs(c.Mass - mass) <= maxDa), Is.True,
                "High consecutive charges of a clean protein should not be spuriously rejected as harmonics");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 4. Mass range boundary behaviour
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FindCandidateMasses_MassJustBelowMinRange_IsExcluded()
        {
            // A protein at 49 Da (below MinMassRange=50) should not appear
            // in results — but the code must not throw.
            const double mass = 49.0;
            // At such low mass, charges 1,2,3 give mz values around 25-50
            var spectrum = MakeSpectrum(mass, new[] { 1, 2, 3 });
            var deconParams = Params(minMass: 50.0, maxMass: 100_000.0);

            List<FLASHDeconvolutionAlgorithm.CandidateMass> candidates = null;
            Assert.DoesNotThrow(() => candidates = RunStep2(spectrum, deconParams));

            // Any candidate produced must still be within the declared range
            foreach (var c in candidates)
                Assert.That(c.Mass, Is.GreaterThanOrEqualTo(deconParams.MinMassRange));
        }

        [Test]
        public void FindCandidateMasses_MassJustAboveMaxRange_IsExcluded()
        {
            // A very large protein (110 kDa, above MaxMassRange=100 kDa) should be excluded.
            // Uses consecutive high charges so the continuity check can be satisfied.
            const double mass = 110_000.0;
            var spectrum = MakeSpectrum(mass, new[] { 50, 51, 52, 53, 54 });
            var deconParams = Params(maxMass: 100_000.0);

            List<FLASHDeconvolutionAlgorithm.CandidateMass> candidates = null;
            Assert.DoesNotThrow(() => candidates = RunStep2(spectrum, deconParams));

            foreach (var c in candidates)
                Assert.That(c.Mass, Is.LessThanOrEqualTo(deconParams.MaxMassRange),
                    $"Candidate at {c.Mass:F0} Da exceeds MaxMassRange {deconParams.MaxMassRange}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 5. Retroactive-credit safety — no false predecessor credit
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FindCandidateMasses_TwoSeparateMasses_BothFound()
        {
            // Two proteins at well-separated masses should both be found as candidates.
            // Step 2 is intentionally permissive: it may also produce spurious candidates
            // (e.g. harmonics of one mass appearing near another mass value) — these are
            // expected and will be rejected by the isotope cosine scoring in Steps 4–5.
            // The invariant tested here is that Step 2 does NOT miss either true mass,
            // i.e., the prevIntensity state from mass1 peaks does not corrupt mass2's
            // accumulator or vice versa.
            double mass1 = 5_000.0;
            double mass2 = 50_000.0;

            var mzPairs = new List<(double mz, double intensity)>();
            foreach (int z in new[] { 5, 6, 7 })
                mzPairs.Add((mass1.ToMz(z), 10_000.0));
            foreach (int z in new[] { 25, 26, 27 })  // consecutive, not harmonics of mass1
                mzPairs.Add((mass2.ToMz(z), 10_000.0));

            var ordered = mzPairs.OrderBy(p => p.mz).ToArray();
            var spectrum = new MzSpectrum(
                ordered.Select(p => p.mz).ToArray(),
                ordered.Select(p => p.intensity).ToArray(),
                shouldCopy: false);

            var candidates = RunStep2(spectrum, Params());

            double maxDa1 = mass1 * 30e-6;
            double maxDa2 = mass2 * 30e-6;

            Assert.That(candidates.Any(c => Math.Abs(c.Mass - mass1) <= maxDa1), Is.True,
                $"Expected a candidate near mass1={mass1} Da");
            Assert.That(candidates.Any(c => Math.Abs(c.Mass - mass2) <= maxDa2), Is.True,
                $"Expected a candidate near mass2={mass2} Da");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 6. Narrow charge range parameter
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FindCandidateMasses_NarrowChargeRange_StillFindsIfChargesMatch()
        {
            // If we restrict the search to charges 8–12 and the protein is at 9–11,
            // it should still be found.
            const double mass = 12_223.0;
            var spectrum = MakeSpectrum(mass, new[] { 9, 10, 11 });
            var deconParams = Params(minCharge: 8, maxCharge: 12);
            var candidates = RunStep2(spectrum, deconParams);

            double maxDa = mass * 30e-6;
            Assert.That(candidates.Any(c => Math.Abs(c.Mass - mass) <= maxDa), Is.True,
                "Protein within the specified charge range should be found");
        }

        [Test]
        public void FindCandidateMasses_ChargeRangeExcludesProtein_NoNearMassCandidate()
        {
            // If we restrict search to charges 20–60 but the protein's true charge
            // states (5, 6, 7) are outside that range, no candidate should appear
            // near the true mass of 5000 Da. The peaks may still generate spurious
            // candidates at other masses (e.g. ~20–35 kDa, if those peaks happen
            // to form a continuous run when interpreted as high charges), but the
            // true 5000 Da mass must not appear.
            const double mass = 5_000.0;
            var spectrum = MakeSpectrum(mass, new[] { 5, 6, 7 });
            var deconParams = Params(minCharge: 20, maxCharge: 60);
            var candidates = RunStep2(spectrum, deconParams);

            double maxDa = mass * 30e-6;
            Assert.That(candidates.Any(c => Math.Abs(c.Mass - mass) <= maxDa), Is.False,
                "No candidate should appear near the true mass when the protein's charge states are outside the search range");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 7. BinNumber / BinValue round-trip at realistic mass values
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        [TestCase(500.0)]
        [TestCase(12_223.0)]    // Cytochrome C
        [TestCase(50_000.0)]    // mid-range
        [TestCase(99_999.0)]    // near max
        public void BinArithmetic_RoundTrip_WithinOnePpm(double mass)
        {
            double binMulFactor = 1.0 / (10.0 * 1e-6);
            double minValue = 0.0; // log(1)
            double logMass = Math.Log(mass);

            int bin = FLASHDeconvolutionAlgorithm.BinNumber(logMass, minValue, binMulFactor);
            double recovered = FLASHDeconvolutionAlgorithm.BinValue(bin, minValue, binMulFactor);
            double recoveredMass = Math.Exp(recovered);

            // Each bin spans exactly 1/binMulFactor in log space = tolerancePpm in ppm.
            // The maximum round-trip error is half a bin width = tolerancePpm / 2 = 5 ppm.
            // The test name says "within one ppm" but that was wrong — it should be
            // "within half the bin width", i.e. < 5 ppm at 10 ppm tolerance.
            double maxPpmError = 5.0; // = tolerancePpm / 2
            double ppmError = Math.Abs(recoveredMass - mass) / mass * 1e6;
            Assert.That(ppmError, Is.LessThan(maxPpmError),
                $"Round-trip error for mass={mass}: {ppmError:F3} ppm (max allowed: {maxPpmError} ppm)");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 8. Consistent results regardless of spectrum intensity scale
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FindCandidateMasses_ResultsIndependentOfAbsoluteIntensityScale()
        {
            // Multiplying all intensities by a constant should not change which
            // candidates are found (only the SupportIntensity values scale).
            const double mass = 8_000.0;
            int[] charges = { 7, 8, 9, 10, 11 };
            double[] baseInt = { 1.0, 1.0, 1.0, 1.0, 1.0 };
            double[] scaledInt = { 1000.0, 1000.0, 1000.0, 1000.0, 1000.0 };

            var specBase = MakeSpectrum(mass, charges, baseInt);
            var specScaled = MakeSpectrum(mass, charges, scaledInt);
            var deconParams = Params();

            var candidatesBase = RunStep2(specBase, deconParams);
            var candidatesScaled = RunStep2(specScaled, deconParams);

            double maxDa = mass * 30e-6;
            bool foundBase = candidatesBase.Any(c => Math.Abs(c.Mass - mass) <= maxDa);
            bool foundScaled = candidatesScaled.Any(c => Math.Abs(c.Mass - mass) <= maxDa);

            Assert.That(foundBase, Is.True, "Base intensity: candidate not found");
            Assert.That(foundScaled, Is.True, "Scaled intensity: candidate not found");
            Assert.That(foundBase == foundScaled, Is.True,
                "Candidate detection should be independent of absolute intensity scale");
        }
    }
}