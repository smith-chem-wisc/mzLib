using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using static Test.MassSpectrometryTests.Deconvolution.DeconvolutionTestHelpers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Unit tests for the spectrum-aware path of <see cref="DeconvolutionScorer"/>:
    /// the <see cref="DeconvolutionScorer.ComputeFeatures(IsotopicEnvelope, AverageResidue, MzSpectrum)"/>
    /// overload, the two new features (LocalSignalToNoise, CompetingPeakRatio),
    /// <see cref="DeconvolutionScorer.ComputeScoreWithSpectrumContext"/>, and the new
    /// three-arg <see cref="DeconvolutionScorer.ScoreEnvelope(IsotopicEnvelope, AverageResidue, MzSpectrum)"/>
    /// overload.
    ///
    /// Test envelopes are built via the shared <see cref="DeconvolutionTestHelpers.BuildPerfectEnvelope()"/>
    /// helper so synthetic-envelope construction matches the envelope-only test fixture
    /// (<c>TestDeconvolutionScorerUnit</c>) exactly. Spectra are built locally via
    /// <see cref="BuildSpectrum"/>.
    /// </summary>
    [TestFixture]
    public sealed class TestDeconvolutionScorerSpectrumAware
    {
        // ── Local helpers ─────────────────────────────────────────────────────

        /// <summary>
        /// Builds an <see cref="MzSpectrum"/> from a list of (mz, intensity) peaks.
        /// Sorts by m/z ascending (the spectrum class assumes XArray is sorted).
        /// </summary>
        private static MzSpectrum BuildSpectrum(IEnumerable<(double mz, double intensity)> peaks)
        {
            var sorted = peaks.OrderBy(p => p.mz).ToArray();
            double[] mz = sorted.Select(p => p.mz).ToArray();
            double[] intensity = sorted.Select(p => p.intensity).ToArray();
            return new MzSpectrum(mz, intensity, false);
        }

        /// <summary>
        /// Builds a clean spectrum: every peak from the envelope, plus optional extra peaks
        /// (noise, competing peaks, etc.).
        /// </summary>
        private static MzSpectrum SpectrumFromEnvelope(IsotopicEnvelope env,
            IEnumerable<(double mz, double intensity)> extraPeaks = null)
        {
            var all = env.Peaks.ToList();
            if (extraPeaks != null) all.AddRange(extraPeaks);
            return BuildSpectrum(all);
        }

        // ── Argument validation ───────────────────────────────────────────────

        [Test]
        public void ComputeFeatures_SpectrumOverload_DefaultsForNullSpectrumThrows()
        {
            var env = BuildPerfectEnvelope();
            var spec = SpectrumFromEnvelope(env);

            var exEnv = Assert.Throws<ArgumentNullException>(
                () => DeconvolutionScorer.ComputeFeatures(null, Model, spec));
            Assert.That(exEnv.ParamName, Is.EqualTo("envelope"));

            var exModel = Assert.Throws<ArgumentNullException>(
                () => DeconvolutionScorer.ComputeFeatures(env, null, spec));
            Assert.That(exModel.ParamName, Is.EqualTo("model"));

            var exSpec = Assert.Throws<ArgumentNullException>(
                () => DeconvolutionScorer.ComputeFeatures(env, Model, (MzSpectrum)null));
            Assert.That(exSpec.ParamName, Is.EqualTo("spectrum"));
        }

        // ── Envelope-only feature preservation ────────────────────────────────

        [Test]
        public void ComputeFeatures_SpectrumOverload_PreservesEnvelopeOnlyFeatureValues()
        {
            // Perfect envelope embedded in a clean spectrum (envelope's own peaks only).
            var env = BuildPerfectEnvelope();
            var spec = SpectrumFromEnvelope(env);

            var envOnly = DeconvolutionScorer.ComputeFeatures(env, Model);
            var withSpec = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            // The four envelope-only fields must be byte-equal — the spectrum-aware overload
            // is required by contract to forward them verbatim.
            Assert.That(withSpec.AveragineCosineSimilarity,
                Is.EqualTo(envOnly.AveragineCosineSimilarity).Within(1e-12));
            Assert.That(withSpec.AvgPpmError,
                Is.EqualTo(envOnly.AvgPpmError).Within(1e-12));
            Assert.That(withSpec.PeakCompleteness,
                Is.EqualTo(envOnly.PeakCompleteness).Within(1e-12));
            Assert.That(withSpec.IntensityRatioConsistency,
                Is.EqualTo(envOnly.IntensityRatioConsistency).Within(1e-12));
        }

        [Test]
        public void ComputeFeatures_SpectrumOverload_HasSpectrumFeaturesIsTrue()
        {
            var env = BuildPerfectEnvelope();
            var spec = SpectrumFromEnvelope(env);

            var f = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            Assert.That(f.HasSpectrumFeatures, Is.True,
                "Spectrum-aware ComputeFeatures must produce HasSpectrumFeatures == true");
            Assert.That(double.IsFinite(f.LocalSignalToNoise), Is.True,
                $"LocalSignalToNoise must be finite, got {f.LocalSignalToNoise}");
            Assert.That(double.IsFinite(f.CompetingPeakRatio), Is.True,
                $"CompetingPeakRatio must be finite, got {f.CompetingPeakRatio}");
        }

        [Test]
        public void ComputeFeatures_NoSpectrum_HasSpectrumFeaturesIsFalse()
        {
            var env = BuildPerfectEnvelope();

            var f = DeconvolutionScorer.ComputeFeatures(env, Model);

            Assert.That(f.HasSpectrumFeatures, Is.False,
                "Envelope-only ComputeFeatures must produce HasSpectrumFeatures == false");
            Assert.That(double.IsNaN(f.LocalSignalToNoise), Is.True,
                "LocalSignalToNoise must be NaN under the envelope-only path");
            Assert.That(double.IsNaN(f.CompetingPeakRatio), Is.True,
                "CompetingPeakRatio must be NaN under the envelope-only path");
        }

        // ── LocalSignalToNoise behaviour ──────────────────────────────────────

        [Test]
        public void LocalSnr_CleanSpectrum_HighSnr()
        {
            // Perfect envelope. Surround each envelope peak with non-envelope peaks far below
            // the envelope's mean peak intensity so the local noise floor is well below every
            // envelope peak. We deliberately use a small fraction (0.1%) rather than a literal
            // 1% so the > 50 bound is robust to whatever specific intensity distribution
            // BuildPerfectEnvelope happens to produce.
            var env = BuildPerfectEnvelope();
            double meanEnvIntensity = env.Peaks.Average(p => p.intensity);
            double noiseLevel = 0.001 * meanEnvIntensity;

            var noisePeaks = new List<(double mz, double intensity)>();
            foreach (var (mz, _) in env.Peaks)
            {
                // Place noise peaks 0.7 and 1.3 Da off — inside the SNR window (±2.0 Da)
                // and outside the competing-peak window (±0.5 Da) and the 10 ppm match
                // tolerance. Picked at non-integer multiples of typical isotope spacing
                // so they don't coincidentally collide with neighbouring envelope peaks.
                noisePeaks.Add((mz - 1.3, noiseLevel));
                noisePeaks.Add((mz - 0.7, noiseLevel));
                noisePeaks.Add((mz + 0.7, noiseLevel));
                noisePeaks.Add((mz + 1.3, noiseLevel));
            }
            var spec = SpectrumFromEnvelope(env, noisePeaks);

            var f = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            Assert.That(f.LocalSignalToNoise, Is.GreaterThan(50.0),
                $"Clean spectrum should yield SNR > 50. Got {f.LocalSignalToNoise:F2}");
        }

        [Test]
        public void LocalSnr_NoisySpectrum_LowSnr()
        {
            // Same envelope. Pad the spectrum with non-envelope peaks at intensity ≈ the
            // envelope's mean peak intensity. The median per-peak SNR collapses toward 1
            // and should be below 5.
            var env = BuildPerfectEnvelope();
            double meanEnvIntensity = env.Peaks.Average(p => p.intensity);
            double noiseLevel = meanEnvIntensity;

            var noisePeaks = new List<(double mz, double intensity)>();
            foreach (var (mz, _) in env.Peaks)
            {
                noisePeaks.Add((mz - 1.3, noiseLevel));
                noisePeaks.Add((mz - 0.7, noiseLevel));
                noisePeaks.Add((mz + 0.7, noiseLevel));
                noisePeaks.Add((mz + 1.3, noiseLevel));
            }
            var spec = SpectrumFromEnvelope(env, noisePeaks);

            var f = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            Assert.That(f.LocalSignalToNoise, Is.LessThan(5.0),
                $"Noisy spectrum should yield SNR < 5. Got {f.LocalSignalToNoise:F2}");
        }

        // ── CompetingPeakRatio behaviour ──────────────────────────────────────

        [Test]
        public void CompetingPeakRatio_NoCompetingPeaks_RatioIsOne()
        {
            // Use charge=1 so envelope peaks are spaced ~1 Da apart (> 0.5 Da CompetingPeakWindowDa).
            // This means each matched position's window contains only its own envelope peak —
            // none of its neighbors. With the spectrum containing only envelope peaks, the
            // apex search excludes them all → apexInWindow = 0 → every matched position is
            // trivially the local apex.
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 1);
            var spec = SpectrumFromEnvelope(env);

            var f = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            Assert.That(f.CompetingPeakRatio, Is.EqualTo(1.0).Within(1e-12),
                $"With no competing peaks, ratio should be 1.0. Got {f.CompetingPeakRatio:F6}");
        }

        [Test]
        public void CompetingPeakRatio_CompetingPeakAtEveryPosition_RatioIsZero()
        {
            // charge=1 so each matched position's window only contains its own envelope peak.
            // Place a competing peak 0.1 Da from each envelope peak at 10× intensity. 0.1 Da
            // is inside the 0.5 Da competing window and far outside the 10 ppm matching
            // tolerance, so each is a non-envelope competitor visible to apex search.
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 1);
            var competing = env.Peaks.Select(p => (p.mz + 0.1, p.intensity * 10.0)).ToList();
            var spec = SpectrumFromEnvelope(env, competing);

            var f = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            Assert.That(f.CompetingPeakRatio, Is.EqualTo(0.0).Within(1e-12),
                $"With a taller competing peak at every matched position, ratio should be 0.0. Got {f.CompetingPeakRatio:F6}");
        }

        [Test]
        public void CompetingPeakRatio_CompetingAtHalfThePositions_RatioNearHalf()
        {
            // charge=1 again. Adding a competing peak only at every other env peak's m/z + 0.1 Da
            // means roughly half the matched positions have a competing peak in their window
            // (the half whose env-peak index is even) and the other half have none. The exact
            // ratio depends on parity alignment between env.Peaks indices and the 1%+ matched
            // subset, but it lands inside [0.4, 0.6] for any reasonable Averagine envelope.
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 1);
            var competing = env.Peaks
                .Select((p, i) => (p, i))
                .Where(x => x.i % 2 == 0)
                .Select(x => (x.p.mz + 0.1, x.p.intensity * 10.0))
                .ToList();
            var spec = SpectrumFromEnvelope(env, competing);

            var f = DeconvolutionScorer.ComputeFeatures(env, Model, spec);

            Assert.That(f.CompetingPeakRatio, Is.InRange(0.4, 0.6),
                $"With competing peaks at ~half the matched positions, ratio should be near 0.5. Got {f.CompetingPeakRatio:F4}");
        }

        // ── ComputeScoreWithSpectrumContext ───────────────────────────────────

        [Test]
        public void ComputeScoreWithSpectrumContext_RequiresSpectrumFeatures()
        {
            // Four-arg constructor → spectrum-aware fields are NaN → must throw.
            var f = new EnvelopeScoreFeatures(0.95, 1.5, 0.95, 0.95);

            var ex = Assert.Throws<ArgumentException>(
                () => DeconvolutionScorer.ComputeScoreWithSpectrumContext(f));
            Assert.That(ex.ParamName, Is.EqualTo("features"));
        }

        [Test]
        public void ComputeScoreWithSpectrumContext_HighQuality_AbovePointFive()
        {
            var f = new EnvelopeScoreFeatures(
                averagineCosineSimilarity: 0.95,
                avgPpmError: 1.5,
                peakCompleteness: 0.95,
                intensityRatioConsistency: 0.95,
                localSignalToNoise: 20.0,
                competingPeakRatio: 1.0);

            double score = DeconvolutionScorer.ComputeScoreWithSpectrumContext(f);

            Assert.That(score, Is.GreaterThan(0.5),
                $"High-quality spectrum-aware features should score > 0.5. Got {score:F4}");
            Assert.That(score, Is.LessThanOrEqualTo(1.0));
        }

        [Test]
        public void ComputeScoreWithSpectrumContext_LowQuality_BelowPointFive()
        {
            var f = new EnvelopeScoreFeatures(
                averagineCosineSimilarity: 0.2,
                avgPpmError: 25.0,
                peakCompleteness: 0.2,
                intensityRatioConsistency: 0.1,
                localSignalToNoise: 0.5,
                competingPeakRatio: 0.0);

            double score = DeconvolutionScorer.ComputeScoreWithSpectrumContext(f);

            Assert.That(score, Is.LessThan(0.5),
                $"Low-quality spectrum-aware features should score < 0.5. Got {score:F4}");
            Assert.That(score, Is.GreaterThanOrEqualTo(0.0));
        }

        [Test]
        public void ComputeScoreWithSpectrumContext_BoundedZeroOne(
            [Values(0.0, 0.5, 1.0)] double cosine,
            [Values(0.0, 0.5, 1.0)] double completeness,
            [Values(0.0, 0.5, 1.0)] double ratio,
            [Values(0.0, 5.0, 50.0)] double snr,
            [Values(0.0, 0.5, 1.0)] double competing)
        {
            var f = new EnvelopeScoreFeatures(
                averagineCosineSimilarity: cosine,
                avgPpmError: 5.0,
                peakCompleteness: completeness,
                intensityRatioConsistency: ratio,
                localSignalToNoise: snr,
                competingPeakRatio: competing);

            double score = DeconvolutionScorer.ComputeScoreWithSpectrumContext(f);

            Assert.That(score, Is.InRange(0.0, 1.0),
                $"Score must be in [0,1]. cosine={cosine} comp={completeness} ratio={ratio} snr={snr} competing={competing} → {score:F6}");
        }

        // ── ScoreEnvelope(env, model, spectrum) overload ──────────────────────

        [Test]
        public void ScoreEnvelope_WithSpectrum_ReturnsSpectrumAwareScore()
        {
            var env = BuildPerfectEnvelope();
            var spec = SpectrumFromEnvelope(env);

            double convenience = DeconvolutionScorer.ScoreEnvelope(env, Model, spec);
            double explicitChain = DeconvolutionScorer.ComputeScoreWithSpectrumContext(
                DeconvolutionScorer.ComputeFeatures(env, Model, spec));

            Assert.That(convenience, Is.EqualTo(explicitChain).Within(1e-12),
                $"ScoreEnvelope(env, model, spec) must equal ComputeScoreWithSpectrumContext(ComputeFeatures(env, model, spec)). Got {convenience} vs {explicitChain}");
        }

        // ── DeconvoluteWithGenericScoring uses the spectrum-aware path ────────

        [Test]
        public void DeconvoluteWithGenericScoring_PopulatesScoreFromSpectrumAwarePath()
        {
            // Part 1: every yielded envelope has GenericScore.HasValue == true.
            var env = BuildPerfectEnvelope(monoMass: TestMass, charge: 2);
            double[] mzArr = env.Peaks.Select(p => p.mz).ToArray();
            double[] yArr = env.Peaks.Select(p => p.intensity).ToArray();
            var clean = new MzSpectrum(mzArr, yArr, false);

            var deconParams = new ClassicDeconvolutionParameters(
                minCharge: 1, maxCharge: 10, deconPpm: 20, intensityRatio: 3);

            var cleanResults = Deconvoluter.DeconvoluteWithGenericScoring(clean, deconParams).ToList();
            Assert.That(cleanResults, Is.Not.Empty, "Clean spectrum should yield at least one envelope");
            foreach (var e in cleanResults)
            {
                Assert.That(e.HasGenericScore, Is.True,
                    "Every yielded envelope must have GenericScore set");
                Assert.That(e.GenericScore!.Value, Is.InRange(0.0, 1.0));
            }

            // Part 2 (stronger — wire-through proof, harmonics-free).
            //
            // We avoid running Classic on a noisy spectrum: noise peaks placed at non-trivial
            // m/z offsets around envelope peaks tend to land on harmonic isotope spacings at
            // some explored charge state (e.g. ±0.7 Da is a near-integer multiple of the z=10
            // step), and Classic can produce different envelopes from clean vs noisy spectra
            // for that reason — confounding any direct comparison of the yielded GenericScore.
            //
            // Instead, take the envelope Classic actually yielded from the clean spectrum and:
            //   (a) verify its GenericScore equals ScoreEnvelope(env, model, cleanSpec) computed
            //       directly — proving the deconvoluter is wiring this exact spectrum into the
            //       scorer.
            //   (b) verify ScoreEnvelope(env, model, cleanSpec) differs from
            //       ScoreEnvelope(env, model, syntheticNoisySpec) — proving the spectrum-aware
            //       scorer is sensitive to spectral context (and not, for example, ignoring the
            //       spectrum and falling back to envelope-only).
            // Together these prove DeconvoluteWithGenericScoring is on the spectrum-aware path.
            var yielded = cleanResults[0];
            double directCleanScore =
                DeconvolutionScorer.ScoreEnvelope(yielded, deconParams.AverageResidueModel, clean);
            Assert.That(yielded.GenericScore!.Value, Is.EqualTo(directCleanScore).Within(1e-12),
                "Yielded GenericScore must match a direct spectrum-aware ScoreEnvelope call against the same spectrum.");

            // Build a synthetic noisy spectrum (NOT fed to the deconvoluter — only used to score
            // the already-yielded envelope). Place noise inside the SNR window (±2 Da) but
            // outside the competing-peak window (±0.5 Da) and outside the 10 ppm match tolerance.
            // Harmonics inside Classic are irrelevant here because we don't run Classic on this
            // spectrum.
            double meanEnvIntensity = yielded.Peaks.Average(p => p.intensity);
            var noisePeaks = new List<(double mz, double intensity)>();
            foreach (var (mz, _) in yielded.Peaks)
            {
                noisePeaks.Add((mz - 1.3, meanEnvIntensity));
                noisePeaks.Add((mz - 0.7, meanEnvIntensity));
                noisePeaks.Add((mz + 0.7, meanEnvIntensity));
                noisePeaks.Add((mz + 1.3, meanEnvIntensity));
            }
            var combined = yielded.Peaks.Concat(noisePeaks).OrderBy(p => p.mz).ToArray();
            var noisy = new MzSpectrum(
                combined.Select(p => p.mz).ToArray(),
                combined.Select(p => p.intensity).ToArray(),
                false);

            double directNoisyScore =
                DeconvolutionScorer.ScoreEnvelope(yielded, deconParams.AverageResidueModel, noisy);
            Assert.That(Math.Abs(directCleanScore - directNoisyScore), Is.GreaterThan(1e-9),
                $"Spectrum-aware scorer must produce different scores for clean vs noisy spectra of the same envelope. " +
                $"Clean: {directCleanScore:F6}, Noisy: {directNoisyScore:F6}");
        }
    }
}
