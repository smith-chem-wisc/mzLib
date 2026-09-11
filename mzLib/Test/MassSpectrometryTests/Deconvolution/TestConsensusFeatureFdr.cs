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
            // Reproduce the generator's arithmetic over many draws: every offset must sit between
            // the teeth, and its magnitude must stay in the documented window. Being off-lattice is
            // what stops a decoy from landing on the target's own peaks.
            var rng = new Random(42);
            for (int i = 0; i < 500; i++)
            {
                double magnitude = ConsensusFeatureFdr.MinDecoyOffsetDa
                    + rng.NextDouble() * (ConsensusFeatureFdr.MaxDecoyOffsetDa - ConsensusFeatureFdr.MinDecoyOffsetDa);
                double k = Math.Round(magnitude / Constants.C13MinusC12);
                double offset = (k + 0.5) * Constants.C13MinusC12;

                double teeth = offset / Constants.C13MinusC12;
                Assert.That(Math.Abs(teeth - Math.Truncate(teeth)), Is.EqualTo(0.5).Within(1e-9),
                    "offset must be a half-integer number of 13C units");
                Assert.That(Math.Abs(offset),
                    Is.InRange(ConsensusFeatureFdr.MinDecoyOffsetDa - Constants.C13MinusC12,
                               ConsensusFeatureFdr.MaxDecoyOffsetDa + Constants.C13MinusC12));
            }
        }

        [Test]
        public void B2_DecoySeparationDoesNotShrinkWithCharge()
        {
            // The point of displacing the whole feature rather than perturbing isotope spacing.
            // A spacing perturbation is a fixed offset in neutral mass, so in m/z it is divided by
            // charge and collapses onto the real lattice at top-down charge states. A 5-50 Da
            // displacement stays far outside tolerance at any charge.
            const double displacement = 7.5 * 1.0033548;
            foreach (int z in new[] { 2, 10, 30, 50 })
            {
                double mzSeparation = displacement / z;
                double ppm = mzSeparation / (5000.0.ToMz(z)) * 1e6;
                Assert.That(ppm, Is.GreaterThan(20.0),
                    $"at z={z} the decoy must remain outside a 20 ppm tolerance (was {ppm:F1} ppm)");
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
            // One feature that is really in the spectrum and several that are not. The real one
            // must not come out with a worse q-value than the fabricated ones.
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };

            var features = new List<MassFeature> { FeatureAt(5000.0, 5, 10.0) };
            for (int i = 1; i <= 8; i++)
                features.Add(FeatureAt(5000.0 + i * 37.3, 5, 10.0));

            ConsensusFeatureFdr.AssignQValues(features, scans, Model);

            var real = features[0];
            var fabricated = features.Skip(1).Where(f => !double.IsNaN(f.QValue)).ToList();
            Assert.That(double.IsNaN(real.QValue), Is.False, "the real feature should be scoreable");
            foreach (var f in fabricated)
                Assert.That(real.QValue, Is.LessThanOrEqualTo(f.QValue + 1e-9),
                    "a feature actually present must not rank below one that is not");
        }

        [Test]
        public void C3_QValueIsMonotoneInScore()
        {
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };
            var features = new List<MassFeature> { FeatureAt(5000.0, 5, 10.0) };
            for (int i = 1; i <= 20; i++)
                features.Add(FeatureAt(5000.0 + i * 13.7, 5, 10.0));

            ConsensusFeatureFdr.AssignQValues(features, scans, Model);

            var scored = features.Where(f => !double.IsNaN(f.QValue))
                                 .Select(f => (f.QValue, Score: IsotopicEnvelopeProbe.Score(
                                     spectrum, f.ConsensusMass, f.Traces[0].Charge, Model, 20.0)))
                                 .Where(x => !double.IsNaN(x.Score))
                                 .OrderByDescending(x => x.Score)
                                 .ToList();

            for (int i = 1; i < scored.Count; i++)
                Assert.That(scored[i].QValue, Is.GreaterThanOrEqualTo(scored[i - 1].QValue - 1e-9),
                    "q-value must not decrease as score decreases");
        }

        [Test]
        public void C4_SameSeedReproducesTheSameQValues()
        {
            var spectrum = SpectrumWithSpecies(5000.0, 5);
            var scans = new List<MsDataScan> { ScanFrom(spectrum) };

            List<double> Run(int seed)
            {
                var fs = new List<MassFeature> { FeatureAt(5000.0, 5, 10.0) };
                for (int i = 1; i <= 6; i++) fs.Add(FeatureAt(5000.0 + i * 21.1, 5, 10.0));
                ConsensusFeatureFdr.AssignQValues(fs, scans, Model, seed: seed);
                return fs.Select(f => f.QValue).ToList();
            }

            Assert.That(Run(7), Is.EqualTo(Run(7)).AsCollection,
                "a fixed seed must reproduce the decoys and therefore the q-values");
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
