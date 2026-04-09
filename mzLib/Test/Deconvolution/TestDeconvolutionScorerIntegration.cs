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
    /// Integration tests for the generic deconvolution scoring pipeline:
    /// <see cref="DeconvolutionScorer"/>, <see cref="Deconvoluter.DeconvoluteWithDecoys"/>,
    /// and <see cref="DeconvolutionQValueCalculator"/>.
    ///
    /// Group D — scorer on real algorithm output (Classic only on master; IsoDec requires DLL)
    /// Group E — DeconvoluteWithDecoys produces usable target and decoy distributions
    /// Group F — q-value end-to-end
    /// Group G — DeconvolutionQValueCalculator unit cases (edge conditions)
    /// </summary>
    [TestFixture]
    public sealed class TestDeconvolutionScorerIntegration
    {
        // ── Shared test mass ──────────────────────────────────────────────────
        // ~5 kDa — clean isotope pattern, within Averagine calibration range.
        private const double MonoMass = 5000.0;
        private static readonly int[] Charges = { 8, 9, 10, 11 };
        private const double BaseIntensity = 1e6;
        private static readonly AverageResidue Model = new Averagine();

        // ── Synthetic spectrum helper ─────────────────────────────────────────

        /// <summary>
        /// Builds a noiseless synthetic spectrum from the Averagine model.
        /// Intensities are drawn from the Averagine distribution, but peak positions
        /// use idealized uniform C13 spacing (monoMass + n * C13MinusC12) rather than
        /// the actual Averagine mass positions (sorted[n].First). At higher masses
        /// (roughly >10 kDa), Averagine mass spacing diverges from uniform C13 spacing,
        /// which can cause peaks to fall outside the scorer's 10 ppm matching window.
        /// This means the spectrum is not truly "perfect" relative to the Averagine model
        /// at high masses. To fix, replace the m/z calculation with sorted[n].First.ToMz(z).
        /// </summary>
        private static MzSpectrum MakeSyntheticSpectrum(
            double monoMass = MonoMass,
            int[] charges = null,
            double baseIntens = BaseIntensity)
        {
            charges ??= Charges;

            int avgIdx = Model.GetMostIntenseMassIndex(monoMass);

            double[] rawMasses = Model.GetAllTheoreticalMasses(avgIdx);
            double[] rawIntens = Model.GetAllTheoreticalIntensities(avgIdx);

            var sorted = rawMasses.Zip(rawIntens)
                .OrderBy(p => p.First)
                .ToArray();

            var mzList = new List<double>();
            var intList = new List<double>();

            foreach (int z in charges)
            {
                for (int n = 0; n < sorted.Length; n++)
                {
                    double intensity = baseIntens * sorted[n].Second;
                    if (intensity < baseIntens * 0.001) continue;

                    double mz = (monoMass + n * Constants.C13MinusC12).ToMz(z);
                    mzList.Add(mz);
                    intList.Add(intensity);
                }
            }

            var pairs = mzList.Zip(intList).OrderBy(p => p.First).ToArray();
            return new MzSpectrum(
                pairs.Select(p => p.First).ToArray(),
                pairs.Select(p => p.Second).ToArray(),
                shouldCopy: false);
        }

        private static ClassicDeconvolutionParameters DefaultClassic()
            => new ClassicDeconvolutionParameters(
                minCharge: 1,
                maxCharge: 20,
                deconPpm: 10.0,
                intensityRatio: 3.0,
                polarity: Polarity.Positive);

        // ══════════════════════════════════════════════════════════════════════
        // Group D — scorer on real algorithm output
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void D1_Scorer_Classic_EnvelopesHavePositiveScore()
        {
            var spectrum = MakeSyntheticSpectrum();
            var envelopes = Deconvoluter.Deconvolute(spectrum, DefaultClassic()).ToList();

            Assert.That(envelopes, Is.Not.Empty,
                "Classic deconvolution should find at least one envelope on the synthetic spectrum");

            var pairs = DeconvolutionScorer.ScoreEnvelopes(envelopes, Model).ToList();

            foreach (var (env, score) in pairs)
                Assert.That(score, Is.GreaterThan(0.0),
                    $"Generic score must be > 0 for every envelope. Got {score:F4} for mass {env.MonoisotopicMass:F2}");
        }

        [Test]
        public void D2_Scorer_Classic_EnvelopeScoresAreInZeroOneRange()
        {
            var spectrum = MakeSyntheticSpectrum();
            var envelopes = Deconvoluter.Deconvolute(spectrum, DefaultClassic()).ToList();

            foreach (var (_, score) in DeconvolutionScorer.ScoreEnvelopes(envelopes, Model))
                Assert.That(score, Is.InRange(0.0, 1.0),
                    $"Generic score must be in [0, 1]. Got {score:F6}");
        }

        [Test]
        public void D3_Scorer_Classic_BestEnvelopeScoresAbovePointFive()
        {
            var spectrum = MakeSyntheticSpectrum();
            var envelopes = Deconvoluter.Deconvolute(spectrum, DefaultClassic()).ToList();

            Assert.That(envelopes, Is.Not.Empty);

            // The envelope closest to the true mass should score well
            var best = envelopes
                .Select(e => (env: e, ppm: Math.Abs(e.MonoisotopicMass - MonoMass) / MonoMass * 1e6))
                .Where(x => x.ppm < 20.0)
                .OrderBy(x => x.ppm)
                .FirstOrDefault();

            Assert.That(best.env, Is.Not.Null,
                "At least one envelope should be within 20 ppm of the true monoisotopic mass");

            double score = DeconvolutionScorer.ScoreEnvelope(best.env, Model);
            Assert.That(score, Is.GreaterThan(0.5),
                $"The best-matching envelope (ppm error: {best.ppm:F1}) should score > 0.5. Got {score:F4}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // Group E — DeconvoluteWithDecoys
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void E1_DeconvoluteWithDecoys_ReturnsBothLists()
        {
            var spectrum = MakeSyntheticSpectrum();
            var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(spectrum, DefaultClassic());

            Assert.That(targets, Is.Not.Null, "Targets list must not be null");
            Assert.That(decoys, Is.Not.Null, "Decoys list must not be null");
            Assert.That(targets.Count, Is.GreaterThan(0),
                "Target deconvolution should find at least one envelope on the synthetic spectrum");
        }

        [Test]
        public void E2_DeconvoluteWithDecoys_AllScoresInValidRange()
        {
            // Classic deconvolution infers isotope spacing from raw peak-to-peak m/z
            // differences rather than from DecoyIsotopeDistance. The decoy parameter
            // therefore has no effect on Classic output — targets and decoys will be
            // identical envelopes. This test verifies the documented limitation: both
            // lists are non-null, all scores are in [0,1], and no exception is thrown.
            // For algorithms that read DecoyIsotopeDistance (e.g. FLASHDeconv), a
            // separate test can assert that decoy scores are lower on average.
            // See GenericDeconvolutionScorer_SessionWrapUp.md — Known Limitations §4.
            var spectrum = MakeSyntheticSpectrum();
            var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(spectrum, DefaultClassic());

            Assert.That(targets.Count, Is.GreaterThan(0), "Target list must be non-empty");
            Console.WriteLine($"Targets: {targets.Count}  Decoys: {decoys.Count}");

            foreach (var e in targets)
            {
                double score = DeconvolutionScorer.ScoreEnvelope(e, Model);
                Assert.That(score, Is.InRange(0.0, 1.0),
                    $"Target score must be in [0,1]. Got {score:F4}");
            }

            foreach (var e in decoys)
            {
                double score = DeconvolutionScorer.ScoreEnvelope(e, Model);
                Assert.That(score, Is.InRange(0.0, 1.0),
                    $"Decoy score must be in [0,1]. Got {score:F4}");
            }
        }

        [Test]
        public void E3_DeconvoluteWithDecoys_TargetAndDecoyCountsAreReasonable()
        {
            var spectrum = MakeSyntheticSpectrum();
            var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(spectrum, DefaultClassic());

            Assert.That(targets.Count, Is.GreaterThanOrEqualTo(1),
                "Target list must contain at least one envelope");
            Assert.DoesNotThrow(() => _ = decoys.Count,
                "Decoy list must be accessible without throwing");
        }

        [Test]
        public void E4_DeconvoluteWithDecoys_MsDataScanOverload_Works()
        {
            var spectrum = MakeSyntheticSpectrum();
            var scan = new MsDataScan(
                spectrum, 1, 1, true, Polarity.Positive,
                0.0, new MzRange(spectrum.Range.Minimum, spectrum.Range.Maximum),
                null, MZAnalyzerType.Orbitrap, spectrum.SumOfAllY,
                null, null, "scan=1");

            (List<IsotopicEnvelope> targets, List<IsotopicEnvelope> decoys) result = default;
            Assert.DoesNotThrow(
                () => result = Deconvoluter.DeconvoluteWithDecoys(scan, DefaultClassic()),
                "MsDataScan overload of DeconvoluteWithDecoys must not throw");

            Assert.That(result.targets, Is.Not.Null);
            Assert.That(result.decoys, Is.Not.Null);
        }

        // ══════════════════════════════════════════════════════════════════════
        // Group F — Q-value end-to-end
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void F1_QValueCalculation_RoundTrip_ProducesValidQValues()
        {
            var spectrum = MakeSyntheticSpectrum();
            var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(spectrum, DefaultClassic());

            var targetPairs = DeconvolutionScorer.ScoreEnvelopes(targets, Model).ToList();
            var decoyPairs = DeconvolutionScorer.ScoreEnvelopes(decoys, Model).ToList();

            var results = DeconvolutionQValueCalculator
                .AssignQValues(targetPairs, decoyPairs)
                .ToList();

            Assert.That(results.Count, Is.EqualTo(targets.Count),
                "AssignQValues must return one result per target envelope");

            foreach (var (_, score, qValue) in results)
            {
                Assert.That(qValue, Is.InRange(0.0, 1.0),
                    $"Q-value must be in [0, 1]. Got {qValue:F4} for score {score:F4}");
            }

            // Monotonicity: sort by score descending; q-values must be non-decreasing
            var sorted = results.OrderByDescending(r => r.Score).ToList();
            for (int i = 1; i < sorted.Count; i++)
            {
                Assert.That(sorted[i].QValue, Is.GreaterThanOrEqualTo(sorted[i - 1].QValue - 1e-9),
                    $"Q-values must be monotonically non-decreasing as score decreases. " +
                    $"Index {i}: score {sorted[i].Score:F4} has q={sorted[i].QValue:F4} " +
                    $"but previous q={sorted[i - 1].QValue:F4}");
            }
        }

        [Test]
        public void F2_QValueCalculation_QValuesAreInValidRange()
        {
            // Classic decoys are identical to targets (DecoyIsotopeDistance not used),
            // so q-values will be at or near 1.0 — not a meaningful FDR estimate, but
            // the q-value calculation itself must still produce valid output.
            // A test asserting q-value < 0.5 would only be meaningful for algorithms
            // whose decoy pass produces genuinely different envelopes (e.g. FLASHDeconv).
            var spectrum = MakeSyntheticSpectrum();
            var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(spectrum, DefaultClassic());

            double[] targetScores = targets
                .Select(e => DeconvolutionScorer.ScoreEnvelope(e, Model))
                .ToArray();
            double[] decoyScores = decoys
                .Select(e => DeconvolutionScorer.ScoreEnvelope(e, Model))
                .ToArray();

            double[] qValues = DeconvolutionQValueCalculator.AssignQValues(targetScores, decoyScores);

            Assert.That(qValues.Length, Is.EqualTo(targetScores.Length),
                "AssignQValues must return one q-value per target score");

            foreach (double q in qValues)
                Assert.That(q, Is.InRange(0.0, 1.0),
                    $"All q-values must be in [0,1]. Got {q:F4}");

            Console.WriteLine($"Q-value range: [{qValues.Min():F4}, {qValues.Max():F4}]");
        }

        [Test]
        public void F3_GetFdrThreshold_DoesNotThrow()
        {
            // Classic decoys are identical to targets, so the 1% FDR threshold will
            // be 1.0 (no target passes). This tests that GetFdrThreshold handles
            // that case cleanly without throwing. The threshold value assertion
            // (< 1.0) is only meaningful for algorithms with real decoy separation.
            var spectrum = MakeSyntheticSpectrum();
            var (targets, decoys) = Deconvoluter.DeconvoluteWithDecoys(spectrum, DefaultClassic());

            double[] targetScores = targets
                .Select(e => DeconvolutionScorer.ScoreEnvelope(e, Model))
                .ToArray();
            double[] decoyScores = decoys
                .Select(e => DeconvolutionScorer.ScoreEnvelope(e, Model))
                .ToArray();

            double threshold = 0.0;
            Assert.DoesNotThrow(
                () => threshold = DeconvolutionQValueCalculator.GetFdrThreshold(
                    targetScores, decoyScores, maxFdr: 0.01),
                "GetFdrThreshold must not throw even when no target achieves 1% FDR");

            // Result must be in [0, 1] regardless of whether any target passes
            Assert.That(threshold, Is.InRange(0.0, 1.0),
                $"GetFdrThreshold must return a value in [0,1]. Got {threshold:F4}");

            Console.WriteLine($"1% FDR score threshold: {threshold:F4}");
        }

        // ══════════════════════════════════════════════════════════════════════
        // Group G — DeconvolutionQValueCalculator unit cases (edge conditions)
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void G1_AssignQValues_EmptyTargets_ReturnsEmpty()
        {
            var result = DeconvolutionQValueCalculator.AssignQValues(
                new double[0], new[] { 0.5, 0.3 });

            Assert.That(result, Is.Empty);
        }

        [Test]
        public void G2_AssignQValues_EmptyDecoys_AllQValuesAreZero()
        {
            double[] qValues = DeconvolutionQValueCalculator.AssignQValues(
                new[] { 0.9, 0.7, 0.5 },
                new double[0]);

            Assert.That(qValues.Length, Is.EqualTo(3));
            Assert.That(qValues, Is.All.EqualTo(0.0).Within(1e-10),
                "No decoys → FDR = 0 for all targets");
        }

        [Test]
        public void G3_AssignQValues_AllDecoysHigherThanTargets_AllQValuesAreOne()
        {
            // Every decoy outscores every target → FDR = 100% everywhere
            double[] qValues = DeconvolutionQValueCalculator.AssignQValues(
                new[] { 0.1, 0.2, 0.3 },   // targets — all low
                new[] { 0.8, 0.9, 1.0 });  // decoys — all high

            Assert.That(qValues.Length, Is.EqualTo(3));
            foreach (double q in qValues)
                Assert.That(q, Is.EqualTo(1.0).Within(1e-10),
                    "When all decoys outscore all targets, q-values should be 1.0");
        }

        [Test]
        public void G4_AssignQValues_MonotonicallyEnforced()
        {
            // Construct a case where naive ordering would violate monotonicity.
            // Targets: 0.95, 0.80, 0.70, 0.50
            // Decoys:  0.85 (only one decoy, between the first and second target)
            // Raw q at threshold 0.95: 0/1 = 0.0
            // Raw q at threshold 0.80: 1/2 = 0.5
            // Raw q at threshold 0.70: 1/3 = 0.33 — this must be floored to 0.33, not 0.5
            // After monotone sweep from bottom, q[0.95] = min(0.0, ...) = 0.0, q[0.80] = 0.33...
            double[] targets = { 0.95, 0.80, 0.70, 0.50 };
            double[] decoys = { 0.85 };

            double[] qValues = DeconvolutionQValueCalculator.AssignQValues(targets, decoys);

            Assert.That(qValues.Length, Is.EqualTo(4));

            // After sorting by score descending: 0.95, 0.80, 0.70, 0.50
            // q must be non-decreasing (i.e. lowest score has highest q)
            var sorted = targets.Zip(qValues)
                .OrderByDescending(p => p.First)
                .Select(p => p.Second)
                .ToList();

            for (int i = 1; i < sorted.Count; i++)
                Assert.That(sorted[i], Is.GreaterThanOrEqualTo(sorted[i - 1] - 1e-9),
                    $"Q-values must be non-decreasing. At index {i}: {sorted[i]:F4} < {sorted[i - 1]:F4}");
        }

        [Test]
        public void G5_GetFdrThreshold_NoPassingTargets_ReturnsOne()
        {
            // All targets have q-value 1.0 (every decoy outscores them)
            double threshold = DeconvolutionQValueCalculator.GetFdrThreshold(
                new[] { 0.1, 0.2 },
                new[] { 0.9, 0.8 },
                maxFdr: 0.01);

            Assert.That(threshold, Is.EqualTo(1.0),
                "GetFdrThreshold should return 1.0 when no target achieves the desired FDR");
        }
    }
}