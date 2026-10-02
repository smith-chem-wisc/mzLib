using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Linq;

namespace Test.Deconvolution
{
    /// <summary>
    /// Tests for <see cref="ShuffledAveragine"/>, the charge-invariant decoy model.
    ///
    /// Tests are organised into groups:
    ///   A — Construction and argument validation
    ///   B — Positions are untouched; intensities are falsified
    ///   C — Determinism and seed sensitivity
    ///   D — Equality and hashing
    ///   E — The charge-invariance property this model exists for
    /// </summary>
    [TestFixture]
    public class TestShuffledAveragine
    {
        /// <summary>
        /// Mirrors <c>AverageResidue.NumAveraginesToGenerate</c>, which is protected and so not
        /// visible here. E3 fails loudly if the real table ever becomes smaller than this.
        /// </summary>
        private const int TableSize = 1500;

        private static readonly Averagine RealModel = new();
        private static readonly ShuffledAveragine ShuffledModel = new(new Averagine());

        /// <summary>An index whose envelope is long enough to actually permute.</summary>
        private static int FirstNonDegenerateIndex()
        {
            for (int i = 0; i < TableSize; i++)
                if (RealModel.GetAllTheoreticalIntensities(i).Length >= 3)
                    return i;
            Assert.Fail("no averagine entry with >= 3 peaks");
            return -1;
        }

        // ── A: Construction ───────────────────────────────────────────────────

        [Test]
        public void A1_DefaultConstruct_DoesNotThrow()
        {
            Assert.That(() => new ShuffledAveragine(new Averagine()), Throws.Nothing);
        }

        [Test]
        public void A2_NullRealModel_Throws()
        {
            Assert.That(() => new ShuffledAveragine(null!), Throws.ArgumentNullException);
        }

        [Test]
        public void A3_SeedIsRecorded()
        {
            Assert.That(new ShuffledAveragine(new Averagine(), 7).ShuffleSeed, Is.EqualTo(7));
        }

        // ── B: Positions untouched, intensities falsified ─────────────────────

        [Test]
        public void B1_MassesAreIdenticalToTheRealModel()
        {
            // The whole point: peak positions stay on the 13C lattice, so there is no
            // mass offset that can be divided away by charge.
            for (int i = 0; i < TableSize; i += 97)
                Assert.That(ShuffledModel.GetAllTheoreticalMasses(i),
                    Is.EqualTo(RealModel.GetAllTheoreticalMasses(i)).AsCollection,
                    $"masses differ at averagine index {i}");
        }

        [Test]
        public void B2_DiffToMonoisotopicIsDelegated()
        {
            for (int i = 0; i < TableSize; i += 97)
                Assert.That(ShuffledModel.GetDiffToMonoisotopic(i),
                    Is.EqualTo(RealModel.GetDiffToMonoisotopic(i)));
        }

        [Test]
        public void B3_IntensitiesArePermutedNotAltered()
        {
            // A permutation preserves the multiset: same values, same total, different order.
            int i = FirstNonDegenerateIndex();
            double[] real = RealModel.GetAllTheoreticalIntensities(i);
            double[] shuffled = ShuffledModel.GetAllTheoreticalIntensities(i);

            Assert.That(shuffled.Length, Is.EqualTo(real.Length));
            Assert.That(shuffled.Sum(), Is.EqualTo(real.Sum()).Within(1e-12),
                "shuffling must not create or destroy intensity");
            Assert.That(shuffled.OrderBy(x => x).ToArray(),
                Is.EqualTo(real.OrderBy(x => x).ToArray()).AsCollection.Within(1e-12),
                "shuffled intensities must be a permutation of the real ones");
        }

        [Test]
        public void B4_ApexStaysAtIndexZero()
        {
            // Index 0 is the apex by convention across the AverageResidue arrays. Permuting it
            // would change envelope selection rather than just the decoy's shape.
            for (int i = 0; i < TableSize; i += 97)
            {
                double[] real = RealModel.GetAllTheoreticalIntensities(i);
                double[] shuffled = ShuffledModel.GetAllTheoreticalIntensities(i);
                Assert.That(shuffled[0], Is.EqualTo(real[0]).Within(1e-12),
                    $"apex intensity moved at averagine index {i}");
                Assert.That(shuffled[0], Is.EqualTo(shuffled.Max()).Within(1e-12),
                    $"apex is no longer the tallest peak at averagine index {i}");
            }
        }

        [Test]
        public void B5_ShapeActuallyChangesSomewhere()
        {
            // A decoy that happened to reproduce every real envelope would be useless. At least
            // one sufficiently long envelope must differ from the real model.
            bool anyDifferent = false;
            for (int i = 0; i < TableSize && !anyDifferent; i++)
            {
                if (ShuffledModel.IsDegenerate(i))
                    continue;
                double[] real = RealModel.GetAllTheoreticalIntensities(i);
                double[] shuffled = ShuffledModel.GetAllTheoreticalIntensities(i);
                anyDifferent = real.Where((t, j) => Math.Abs(t - shuffled[j]) > 1e-12).Any();
            }
            Assert.That(anyDifferent, Is.True, "no envelope shape was altered");
        }

        [Test]
        public void B6_DegenerateEnvelopesAreReportedAsSuch()
        {
            // Envelopes with fewer than three peaks have nothing to permute below the apex, so
            // they are not decoys. The model must say so rather than pass them off silently.
            for (int i = 0; i < TableSize; i += 97)
            {
                bool degenerate = RealModel.GetAllTheoreticalIntensities(i).Length < 3;
                Assert.That(ShuffledModel.IsDegenerate(i), Is.EqualTo(degenerate),
                    $"IsDegenerate disagrees with envelope length at index {i}");
            }
        }

        // ── C: Determinism ────────────────────────────────────────────────────

        [Test]
        public void C1_SameSeedGivesSameTable()
        {
            var a = new ShuffledAveragine(new Averagine(), 123);
            var b = new ShuffledAveragine(new Averagine(), 123);
            for (int i = 0; i < TableSize; i += 97)
                Assert.That(a.GetAllTheoreticalIntensities(i),
                    Is.EqualTo(b.GetAllTheoreticalIntensities(i)).AsCollection,
                    $"same seed produced different tables at index {i}");
        }

        [Test]
        public void C2_DifferentSeedsGiveDifferentTables()
        {
            var a = new ShuffledAveragine(new Averagine(), 1);
            var b = new ShuffledAveragine(new Averagine(), 2);
            bool anyDifferent = false;
            for (int i = 0; i < TableSize && !anyDifferent; i++)
            {
                if (a.IsDegenerate(i)) continue;
                anyDifferent = a.GetAllTheoreticalIntensities(i)
                    .Where((t, j) => Math.Abs(t - b.GetAllTheoreticalIntensities(i)[j]) > 1e-12)
                    .Any();
            }
            Assert.That(anyDifferent, Is.True, "two seeds produced identical tables");
        }

        // ── D: Equality ───────────────────────────────────────────────────────

        [Test]
        public void D1_SameSeedModelsAreEqual()
        {
            var a = new ShuffledAveragine(new Averagine(), 5);
            var b = new ShuffledAveragine(new Averagine(), 5);
            Assert.That(a, Is.EqualTo(b));
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void D2_DifferentSeedModelsAreNotEqual()
        {
            Assert.That(new ShuffledAveragine(new Averagine(), 5),
                Is.Not.EqualTo(new ShuffledAveragine(new Averagine(), 6)));
        }

        [Test]
        public void D3_NotEqualToTheRealOrShiftedModel()
        {
            Assert.That(ShuffledModel, Is.Not.EqualTo(RealModel));
            Assert.That(ShuffledModel, Is.Not.EqualTo(new DecoyAveragine(new Averagine())));
        }

        // ── E: The property this model exists for ─────────────────────────────

        [Test]
        public void E1_ShuffledDecoyIsChargeInvariant_ShiftedDecoyIsNot()
        {
            // DecoyAveragine displaces the n-th tooth by a fixed offset in NEUTRAL MASS.
            // Deconvolution matches in m/z, so the matcher sees that offset divided by charge:
            // the decoy converges onto the real lattice as charge rises. ShuffledAveragine moves
            // nothing, so its distance from the real model does not depend on charge at all.
            const double mz = 900.0;                       // where both BU and TD ions sit
            double perTooth = Constants.C13MinusC12 - DecoyAveragine.DefaultDecoyIsotopeSpacing;

            double PpmAtCharge(int z) => (perTooth / z) / mz * 1e6;

            double lowCharge = PpmAtCharge(2);
            double highCharge = PpmAtCharge(30);

            Assert.That(lowCharge, Is.GreaterThan(20.0),
                "at z=2 the shifted decoy should sit outside a 20 ppm tolerance");
            Assert.That(highCharge, Is.LessThan(4.0),
                "at z=30 the shifted decoy should have collapsed inside even a 4 ppm tolerance");
            Assert.That(lowCharge / highCharge, Is.EqualTo(15.0).Within(1e-9),
                "the shifted decoy's separation should scale exactly as 1/z");

            // The shuffled model's masses are charge-independent because they are unchanged.
            int i = FirstNonDegenerateIndex();
            Assert.That(ShuffledModel.GetAllTheoreticalMasses(i),
                Is.EqualTo(RealModel.GetAllTheoreticalMasses(i)).AsCollection,
                "shuffled decoy must not move peaks at all");
        }

        [Test]
        public void E3_MirroredTableSizeStillMatchesTheRealTable()
        {
            // TableSize duplicates a protected constant. If the real table ever shrinks, every
            // sampled loop above would silently index out of range -- catch it here instead.
            Assert.That(() => RealModel.GetAllTheoreticalIntensities(TableSize - 1),
                Throws.Nothing, "the real averagine table is smaller than TableSize");
        }

        [Test]
        public void E2_ShapeDivergenceDoesNotDependOnMass()
        {
            // Sanity: the falsification is present across the mass range, not only for small
            // averagines. Top-down envelopes are long, so if anything the shuffle bites harder.
            int checkedEntries = 0, differing = 0;
            for (int i = 0; i < TableSize; i += 251)
            {
                if (ShuffledModel.IsDegenerate(i)) continue;
                checkedEntries++;
                double[] real = RealModel.GetAllTheoreticalIntensities(i);
                double[] shuffled = ShuffledModel.GetAllTheoreticalIntensities(i);
                if (real.Where((t, j) => Math.Abs(t - shuffled[j]) > 1e-12).Any())
                    differing++;
            }
            Assert.That(checkedEntries, Is.GreaterThan(0), "no non-degenerate entries sampled");
            Assert.That(differing, Is.EqualTo(checkedEntries),
                "every sampled non-degenerate envelope should have been reshaped");
        }
    }
}
