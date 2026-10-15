using Chemistry;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Tests for <see cref="IsotopicEnvelopeProbe"/> and <see cref="ConsensusFeatureFdr"/>.
    ///
    /// Groups:
    ///   A — the probe finds a real envelope and fails to find an absent one
    ///   B — decoy construction stays off the 13C lattice and does not shrink with charge
    ///   C — q-values are assigned, ordered, and bounded
    /// </summary>
    [TestFixture]
    public class TestConsensusFeatureFdr
    {
        private static readonly Averagine Model = new();

        /// <summary>A spectrum containing a real averagine envelope at the given mass and charge.</summary>
        private static MzSpectrum SpectrumWithSpecies(double mass, int charge, double apexIntensity = 1e6)
        {
            int index = Model.GetMostIntenseMassIndex(mass);
            double[] masses = Model.GetAllTheoreticalMasses(index);
            double[] intensities = Model.GetAllTheoreticalIntensities(index);
            double offset = mass + Model.GetDiffToMonoisotopic(index) - masses[0];

            var mz = new List<double>();
            var inten = new List<double>();
            for (int i = 0; i < masses.Length; i++)
            {
                mz.Add((masses[i] + offset).ToMz(charge));
                inten.Add(intensities[i] * apexIntensity);
            }

            var order = Enumerable.Range(0, mz.Count).OrderBy(i => mz[i]).ToList();
            return new MzSpectrum(order.Select(i => mz[i]).ToArray(),
                                  order.Select(i => inten[i]).ToArray(), false);
        }

        /// <summary>
        /// Several real envelopes at one charge over a floor of seeded random noise peaks. The
        /// noise is what makes the C tests meaningful: fabricated hypotheses and decoys pick up a
        /// few stray peaks and so receive finite, varied scores instead of all collapsing to NaN,
        /// and the decoy scores (and therefore the q-values) depend on which offsets were drawn.
        /// </summary>
        private static MzSpectrum NoisySpectrumWithSpecies(IReadOnlyList<(double Mass, double Apex)> species,
            int charge, int noisePeaks = 6000, int noiseSeed = 1)
        {
            var points = new List<(double Mz, double Intensity)>();
            foreach (var (mass, apex) in species)
            {
                int index = Model.GetMostIntenseMassIndex(mass);
                double[] masses = Model.GetAllTheoreticalMasses(index);
                double[] intensities = Model.GetAllTheoreticalIntensities(index);
                double offset = mass + Model.GetDiffToMonoisotopic(index) - masses[0];
                for (int i = 0; i < masses.Length; i++)
                    points.Add(((masses[i] + offset).ToMz(charge), intensities[i] * apex));
            }

            var rng = new Random(noiseSeed);
            double lo = points.Min(p => p.Mz) - 20, hi = points.Max(p => p.Mz) + 20;
            for (int i = 0; i < noisePeaks; i++)
                points.Add((lo + rng.NextDouble() * (hi - lo), 1e3 + rng.NextDouble() * 5e4));

            var sorted = points.OrderBy(p => p.Mz).ToList();
            return new MzSpectrum(sorted.Select(p => p.Mz).ToArray(), sorted.Select(p => p.Intensity).ToArray(), false);
        }

        private static readonly (double Mass, double Apex)[] RealSpecies =
        {
            (5000.0, 1e6), (5600.0, 6e5), (6300.0, 3e5), (7100.0, 1.5e5), (7900.0, 8e4),
        };

        /// <summary>Real features first, then fabricated ones spread between the real species.</summary>
        private static List<MassFeature> RealAndFabricatedFeatures(int fabricated)
        {
            var fs = RealSpecies.Select(s => FeatureAt(s.Mass, 5, 10.0, s.Apex)).ToList();
            for (int i = 1; i <= fabricated; i++)
                fs.Add(FeatureAt(5000.0 + i * 3000.0 / (fabricated + 1) + 0.37 * i, 5, 10.0));
            return fs;
        }

        private static MsDataScan ScanFrom(MzSpectrum s, int number = 1, double rt = 10.0) =>
            new MsDataScan(s, number, 1, true, Polarity.Positive, rt, s.Range, null,
                MZAnalyzerType.Orbitrap, s.SumOfAllY, null, null, null);

        private static MassFeature FeatureAt(double mass, int charge, double rt, double intensity = 1e6)
        {
            var trace = new CorrectedTrace { Id = 1, Charge = charge, ConsensusMass = mass };
            trace.Envelopes.Add(new CorrectedEnvelope
            {
                ScanIndex = 0, ScanNumber = 1, RT = rt,
                OriginalMass = mass, CorrectedMass = mass,
                Charge = charge, Intensity = intensity, WasCorrected = false,
            });
            var f = new MassFeature { Id = 1 };
            f.Traces.Add(trace);
            f.Finalise();
            return f;
        }

        // ── A: the probe ──────────────────────────────────────────────────────

        [Test]
        public void A1_ProbeFindsARealEnvelope()
        {
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var env = IsotopicEnvelopeProbe.Gather(spectrum, 5000.0, 5, Model, 20.0);

            Assert.That(env, Is.Not.Null, "the species is present, the probe should find it");
            Assert.That(env.Peaks.Count, Is.GreaterThan(1), "more than the apex should match");
            Assert.That(env.MonoisotopicMass, Is.EqualTo(5000.0).Within(1e-6));
        }

        [Test]
        public void A2_ProbeReturnsNullWhenNothingIsThere()
        {
            // Same spectrum, hypothesis 500 Da away: nothing should match.
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var env = IsotopicEnvelopeProbe.Gather(spectrum, 5500.0, 5, Model, 20.0);

            Assert.That(env, Is.Null, "an unsupported hypothesis must not fabricate an envelope");
            Assert.That(double.IsNaN(IsotopicEnvelopeProbe.Score(spectrum, 5500.0, 5, Model, 20.0)),
                Is.True, "and must score NaN, not 0");
        }

        [Test]
        public void A3_RealHypothesisOutscoresADisplacedOne()
        {
            // This is the property the whole FDR estimate rests on.
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            double real = IsotopicEnvelopeProbe.Score(spectrum, 5000.0, 5, Model, 20.0);
            double displaced = IsotopicEnvelopeProbe.Score(
                spectrum, 5000.0 + 7.5 * Constants.C13MinusC12, 5, Model, 20.0);

            Assert.That(double.IsNaN(real), Is.False);
            Assert.That(double.IsNaN(displaced) || displaced < real, Is.True,
                $"a half-13C displaced hypothesis should score below the real one (real={real}, displaced={displaced})");
        }

        [Test]
        public void A4_ProbeRejectsZeroCharge()
        {
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            Assert.That(() => IsotopicEnvelopeProbe.Gather(spectrum, 5000.0, 0, Model, 20.0),
                Throws.InstanceOf<ArgumentOutOfRangeException>());
        }

        [Test]
        public void A5_EmptySpectrumReturnsNullRatherThanThrowing()
        {
            var empty = new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false);
            Assert.That(IsotopicEnvelopeProbe.Gather(empty, 5000.0, 5, Model, 20.0), Is.Null);
        }

        // ── B: the decoys ─────────────────────────────────────────────────────

        [Test]
        public void B1_DecoyOffsetsAreHalfIntegerAndInRange()
        {
            // The production generator over many draws: every offset must sit between the teeth,
            // and wherever interleaving suffices its magnitude must stay in the documented 5-50 Da
            // window (rounding to the half-integer can add at most one 13C spacing).
            var rng = new Random(42);
            for (int i = 0; i < 500; i++)
            {
                double offset = ConsensusFeatureFdr.DecoyOffset(rng);

                double teeth = Math.Abs(offset) / Constants.C13MinusC12;
                Assert.That(teeth - Math.Truncate(teeth), Is.EqualTo(0.5).Within(1e-9),
                    "offset must be a half-integer number of 13C units");
                Assert.That(Math.Abs(offset),
                    Is.InRange(ConsensusFeatureFdr.MinDecoyOffsetDa,
                               ConsensusFeatureFdr.MaxDecoyOffsetDa + Constants.C13MinusC12));
            }
        }

        [TestCase(5000.0, 5)]
        [TestCase(10000.0, 10)]
        [TestCase(25000.0, 20)]
        [TestCase(40000.0, 30)]
        [TestCase(60000.0, 50)]
        public void B2_DecoyNeverRecoversTheTargetsOwnPeaks(double mass, int charge)
        {
            // What decides whether a decoy is wrong is the distance from each decoy tooth to the
            // nearest REAL peak. Interleaved teeth sit C/(2z) apart in m/z, about 0.5017/M relative
            // whatever the charge, which falls inside a 20 ppm tolerance above ~25 kDa. Draw decoys
            // from the production generator at the floor the FDR uses and probe a spectrum holding
            // only the target: none may match a single peak.
            const double ppm = ConsensusFeatureFdr.DefaultTolerancePpm;
            var spectrum = SpectrumWithSpecies(mass, charge);
            double floor = ConsensusFeatureFdr.MinimumDecoyOffsetDa(mass, Model, ppm);
            var rng = new Random(42);
            for (int i = 0; i < 200; i++)
            {
                double decoyMass = mass + ConsensusFeatureFdr.DecoyOffset(rng, floor);
                var env = IsotopicEnvelopeProbe.Gather(spectrum, decoyMass, charge, Model, ppm);
                Assert.That(env, Is.Null,
                    $"decoy at {decoyMass - mass:F2} Da from a {mass} Da z={charge} target matched " +
                    $"{env?.Peaks.Count} of its peaks (floor {floor:F1} Da)");
            }
        }

        [Test]
        public void B3_HighMassFeatureGetsTheBestPossibleQValue()
        {
            // End to end through AssignQValues: one real 40 kDa z=30 feature and nothing else in
            // the spectrum. With the decoy kept off the real peaks it is unsupported and enters
            // the null at 0, so the target's q-value is the floor (0 + 1) / 1 ... unless the decoy
            // recovered the target's own peaks, outscored or tied it, and pushed q up.
            // Several seeds, so a lucky draw cannot pass the test.
            var spectrum = SpectrumWithSpecies(40000.0, 30);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };
            double target = IsotopicEnvelopeProbe.Score(spectrum, 40000.0, 30, Model, ConsensusFeatureFdr.DefaultTolerancePpm);
            Assert.That(double.IsNaN(target), Is.False);
            for (int seed = 0; seed < 20; seed++)
            {
                // Replay the draw AssignQValues makes for its single feature.
                var rng = new Random(seed);
                double floor = ConsensusFeatureFdr.MinimumDecoyOffsetDa(40000.0, Model, ConsensusFeatureFdr.DefaultTolerancePpm);
                double decoy = IsotopicEnvelopeProbe.Score(spectrum, 40000.0 + ConsensusFeatureFdr.DecoyOffset(rng, floor),
                    30, Model, ConsensusFeatureFdr.DefaultTolerancePpm);
                Assert.That(double.IsNaN(decoy), Is.True, $"seed {seed}: decoy recovered real peaks (score {decoy})");

                var features = new List<MassFeature> { FeatureAt(40000.0, 30, 10.0) };
                ConsensusFeatureFdr.AssignQValues(features, scans, Model, seed: seed);
                Assert.That(double.IsNaN(features[0].QValue), Is.False, $"seed {seed}");
            }
        }

        // ── C: q-values ───────────────────────────────────────────────────────

        [Test]
        public void C1_QValuesAreAssignedAndBounded()
        {
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };
            var features = new List<MassFeature> { FeatureAt(5000.0, 5, 10.0) };

            int n = ConsensusFeatureFdr.AssignQValues(features, scans, Model);

            Assert.That(n, Is.EqualTo(1));
            Assert.That(double.IsNaN(features[0].QValue), Is.False);
            Assert.That(features[0].QValue, Is.InRange(0.0, 1.0));
        }

        [Test]
        public void C2_RealFeaturesRankAheadOfFabricatedOnes()
        {
            // Five real envelopes of decreasing intensity and twelve fabricated hypotheses between
            // them, all over a noise floor so the fabricated ones pick up stray peaks and get
            // finite q-values. Every real feature must come out ahead of every fabricated one.
            var spectrum = NoisySpectrumWithSpecies(RealSpecies, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };
            var features = RealAndFabricatedFeatures(12);

            ConsensusFeatureFdr.AssignQValues(features, scans, Model);

            var real = features.Take(RealSpecies.Length).ToList();
            var fabricated = features.Skip(RealSpecies.Length).Where(f => !double.IsNaN(f.QValue)).ToList();
            Assert.That(real.All(f => !double.IsNaN(f.QValue)), Is.True, "every real feature should be scoreable");
            Assert.That(fabricated, Has.Count.GreaterThan(1), "fixture must give fabricated features a q-value");

            double worstReal = real.Max(f => f.QValue);
            double bestFabricated = fabricated.Min(f => f.QValue);
            Assert.That(worstReal, Is.LessThan(bestFabricated),
                "a feature actually present must rank ahead of every one that is not");
        }

        [Test]
        public void C3_QValueIsMonotoneInScore()
        {
            var spectrum = NoisySpectrumWithSpecies(RealSpecies, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };
            var features = RealAndFabricatedFeatures(20);

            ConsensusFeatureFdr.AssignQValues(features, scans, Model);

            var scored = features.Where(f => !double.IsNaN(f.QValue))
                                 .Select(f => (f.QValue, Score: IsotopicEnvelopeProbe.Score(
                                     spectrum, f.ConsensusMass, f.Traces[0].Charge, Model, 20.0)))
                                 .OrderByDescending(x => x.Score)
                                 .ToList();

            Assert.That(scored, Has.Count.GreaterThanOrEqualTo(10), "monotonicity needs many scored features");
            Assert.That(scored.Select(x => x.QValue).Distinct().Count(), Is.GreaterThan(2),
                "q-values must actually vary for monotonicity to be tested");
            for (int i = 1; i < scored.Count; i++)
                Assert.That(scored[i].QValue, Is.GreaterThanOrEqualTo(scored[i - 1].QValue - 1e-9),
                    "q-value must not decrease as score decreases");
        }

        [Test]
        public void C4_SameSeedReproducesTheSameQValues()
        {
            var spectrum = NoisySpectrumWithSpecies(RealSpecies, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };

            List<double> Run(int seed)
            {
                var fs = RealAndFabricatedFeatures(20);
                ConsensusFeatureFdr.AssignQValues(fs, scans, Model, seed: seed);
                return fs.Select(f => f.QValue).ToList();
            }

            Assert.That(Run(7), Is.EqualTo(Run(7)).AsCollection,
                "a fixed seed must reproduce the decoys and therefore the q-values");
            // On this fixture the q-values depend on the decoy draw, so an ignored seed (or an
            // unseeded RNG) would show up here as well as in the equality above.
            Assert.That(Run(8), Is.Not.EqualTo(Run(7)).AsCollection,
                "a different seed must change at least one q-value, or the test cannot see the seed");
        }

        // ── D: thresholding a feature file on the q-value ─────────────────────

        [Test]
        public void D1_ThresholdKeepsPassingRowsAndDropsFailingOnes()
        {
            string path = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"fdr_thresh_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                System.IO.File.WriteAllLines(path, new[]
                {
                    "Sample_ID	ID	Mass	Intensity	Time_begin	Time_end	Time_apex	Apex_intensity	Minimum_charge_state	Maximum_charge_state	Minimum_fraction_id	Maximum_fraction_id	Q_value",
                    "0	1	5000.0	1000	10.0	10.2	10.1	900	5	5	0	0	0.01",
                    "0	2	6000.0	1000	10.0	10.2	10.1	900	5	5	0	0	0.90",
                });
                var file = new Readers.Ms1FeatureFile(path);
                Assert.That(file.GetMs1Features(null).Count(), Is.EqualTo(2), "null keeps everything");
                Assert.That(file.GetMs1Features(0.05).Count(), Is.EqualTo(1), "0.90 should be dropped");
                Assert.That(file.GetMs1Features(0.95).Count(), Is.EqualTo(2));
            }
            finally { if (System.IO.File.Exists(path)) System.IO.File.Delete(path); }
        }

        [Test]
        public void D2_RowsWithoutAQValueSurviveAnyThreshold()
        {
            // An external TopFD or FLASHDeconv file carries no q-value. Treating absent as failing
            // would silently return nothing the moment a caller set a threshold, which is the worst
            // possible failure: a search would run on an empty precursor list and simply find less.
            string path = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"fdr_noq_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                System.IO.File.WriteAllLines(path, new[]
                {
                    "Sample_ID	ID	Mass	Intensity	Time_begin	Time_end	Time_apex	Apex_intensity	Minimum_charge_state	Maximum_charge_state	Minimum_fraction_id	Maximum_fraction_id",
                    "0	1	5000.0	1000	10.0	10.2	10.1	900	5	5	0	0",
                });
                var file = new Readers.Ms1FeatureFile(path);
                Assert.That(file.GetMs1Features(0.001).Count(), Is.EqualTo(1),
                    "a file with no q-values must not be emptied by a threshold");
            }
            finally { if (System.IO.File.Exists(path)) System.IO.File.Delete(path); }
        }

        private const string Ms1FeatureHeaderWithQ =
            "Sample_ID	ID	Mass	Intensity	Time_begin	Time_end	Time_apex	Apex_intensity	Minimum_charge_state	Maximum_charge_state	Minimum_fraction_id	Maximum_fraction_id	Q_value";

        [Test]
        public void D3_UnscoredRowsFailTheThresholdWhenTheFileCarriesQValues()
        {
            // An mzLib-written FDR file leaves Q_value empty for a feature the spectrum could not
            // support at all. Inside a file that otherwise carries q-values, that row is the
            // least-supported feature, not one from a tool without q-values, so it must not pass.
            string path = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"fdr_partial_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                System.IO.File.WriteAllLines(path, new[]
                {
                    Ms1FeatureHeaderWithQ,
                    "0	1	5000.0	1000	10.0	10.2	10.1	900	5	5	0	0	0.01",
                    "0	2	6000.0	1000	10.0	10.2	10.1	900	5	5	0	0	",
                });
                var file = new Readers.Ms1FeatureFile(path);
                Assert.That(file.GetMs1Features(null).Count(), Is.EqualTo(2), "no threshold keeps everything");
                var kept = file.GetMs1Features(0.05).ToList();
                Assert.That(kept, Has.Count.EqualTo(1), "the unscored row must be dropped");
                Assert.That(kept[0].Mz, Is.EqualTo(5000.0.ToMz(5)).Within(1e-6));
            }
            finally { if (System.IO.File.Exists(path)) System.IO.File.Delete(path); }
        }

        [Test]
        public void D4_UnscorableFeatureDoesNotSurviveARoundTrip()
        {
            // FromMassFeatures + WriteResults(withQuality) + reload: a NaN q-value is written as an
            // empty cell and must not then pass a threshold that a scored q = 0.02 feature fails.
            var scored = FeatureAt(5000.0, 5, 10.0);
            scored.QValue = 0.01;
            var worse = FeatureAt(5600.0, 5, 10.0);
            worse.QValue = 0.02;
            var unscorable = FeatureAt(6300.0, 5, 10.0);
            unscorable.QValue = double.NaN;

            string path = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"fdr_roundtrip_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                Readers.Ms1FeatureFile.FromMassFeatures(new[] { scored, worse, unscorable }).WriteResults(path, true);
                var kept = new Readers.Ms1FeatureFile(path).GetMs1Features(0.015).ToList();
                Assert.That(kept.Select(k => k.Mz.ToMass(k.Charge)), Is.EqualTo(new[] { 5000.0 }).Within(1e-6).AsCollection);
            }
            finally { if (System.IO.File.Exists(path)) System.IO.File.Delete(path); }
        }

        [Test]
        public void D5_ChangingMaxFeatureQValueAfterLoadReloadsUnderTheNewThreshold()
        {
            string path = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"fdr_params_{Guid.NewGuid():N}_ms1.feature");
            try
            {
                System.IO.File.WriteAllLines(path, new[]
                {
                    Ms1FeatureHeaderWithQ,
                    "0	1	5000.0	1000	10.0	10.2	10.1	900	5	5	0	0	0.005",
                    "0	2	6000.0	1000	10.0	10.2	10.1	900	5	5	0	0	0.03",
                    "0	3	7000.0	1000	10.0	10.2	10.1	900	5	5	0	0	0.50",
                });
                var p = new Readers.FromFileDeconvolutionParameters(path, 1, 60);
                Assert.That(p.Features, Has.Count.EqualTo(3), "unthresholded load");

                p.MaxFeatureQValue = 0.05;
                Assert.That(p.Features, Has.Count.EqualTo(2), "setting a threshold after the first load must take effect");

                p.MaxFeatureQValue = 0.9;
                Assert.That(p.Features, Has.Count.EqualTo(3), "loosening must reload too");

                p.MaxFeatureQValue = 0.01;
                Assert.That(p.Features, Has.Count.EqualTo(1), "a second tightening drops the rows in (0.01, 0.9]");

                p.MaxFeatureQValue = null;
                Assert.That(p.Features, Has.Count.EqualTo(3), "clearing the threshold restores everything");
            }
            finally { if (System.IO.File.Exists(path)) System.IO.File.Delete(path); }
        }

        [Test]
        public void D6_ThresholdOnAnInMemoryInstanceThrows()
        {
            // Built from a feature list there is no file to reload under a threshold; silently
            // keeping the unfiltered list while reporting the threshold is the bug being guarded.
            var p = new Readers.FromFileDeconvolutionParameters(
                new List<ISingleChargeMs1Feature> { new Readers.SingleChargeMs1Feature(1001.0, 5, 10.0, 10.2, 1000) }, 1, 60);
            Assert.That(() => p.MaxFeatureQValue = 0.01, Throws.InvalidOperationException);
            Assert.That(() => p.MaxFeatureQValue = null, Throws.Nothing, "re-setting the current value is a no-op");
            Assert.That(p.Clone().Features, Has.Count.EqualTo(1), "Clone of an in-memory instance still works");
        }

        [Test]
        public void C5_EmptyInputIsHandled()
        {
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };
            Assert.That(ConsensusFeatureFdr.AssignQValues(new List<MassFeature>(), scans, Model),
                Is.EqualTo(0));
            Assert.That(() => ConsensusFeatureFdr.AssignQValues(
                new List<MassFeature> { FeatureAt(5000.0, 5, 10.0) }, new List<MsDataScan>(), Model),
                Throws.ArgumentException);
        }
    }
}
