using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Critical infrastructure tests for the FLASHDeconv deconvolution method.
    /// These tests verify that the enum, parameters class, algorithm skeleton,
    /// and Deconvoluter factory are all correctly wired together.
    ///
    /// Tests do NOT assert on deconvolution quality — the algorithm is a skeleton
    /// that returns empty results. Update assertions in region "Skeleton Behaviour"
    /// once the algorithm is implemented.
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvolutionInfrastructure
    {
        // ── shared spectrum fixtures ──────────────────────────────────────────

        /// <summary>Small positive-mode spectrum with a handful of peaks.</summary>
        private static MzSpectrum SmallPositiveSpectrum() =>
            new MzSpectrum(
                new double[] { 500.0, 500.5, 501.0, 600.0, 601.0 },
                new double[] { 1000.0, 800.0, 400.0, 500.0, 250.0 },
                shouldCopy: false);

        /// <summary>Empty spectrum — zero peaks.</summary>
        private static MzSpectrum EmptySpectrum() =>
            new MzSpectrum(new double[0], new double[0], shouldCopy: false);

        // ══════════════════════════════════════════════════════════════════════
        // 1. Parameters — construction and default values
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FLASHDeconvParameters_DefaultConstructor_DeconvolutionTypeIsCorrect()
        {
            var p = new FLASHDeconvolutionParameters();
            Assert.That(p.DeconvolutionType, Is.EqualTo(DeconvolutionType.FLASHDeconvolution));
        }

        [Test]
        public void FLASHDeconvParameters_DefaultConstructor_ChargeRangeIsCorrect()
        {
            var p = new FLASHDeconvolutionParameters();
            Assert.That(p.MinAssumedChargeState, Is.EqualTo(1));
            Assert.That(p.MaxAssumedChargeState, Is.EqualTo(60));
        }

        [Test]
        public void FLASHDeconvParameters_DefaultConstructor_AllDefaultValuesAreCorrect()
        {
            var p = new FLASHDeconvolutionParameters();

            Assert.That(p.MinIsotopicPeakCount, Is.EqualTo(3));
            Assert.That(p.MaxIsotopicPeakCount, Is.EqualTo(50));
            Assert.That(p.PrecursorHarmonicCount, Is.EqualTo(2));
            Assert.That(p.DeconvolutionTolerancePpm, Is.EqualTo(10.0).Within(1e-10));
            Assert.That(p.MinCosineScore, Is.EqualTo(0.6).Within(1e-10));
            Assert.That(p.IsotopeIntensityRatioThreshold, Is.EqualTo(0.2).Within(1e-10));
            Assert.That(p.MinMassRange, Is.EqualTo(50.0).Within(1e-10));
            Assert.That(p.MaxMassRange, Is.EqualTo(100000.0).Within(1e-10));
        }

        [Test]
        public void FLASHDeconvParameters_DefaultConstructor_PolarityIsPositive()
        {
            var p = new FLASHDeconvolutionParameters();
            Assert.That(p.Polarity, Is.EqualTo(Polarity.Positive));
        }

        [Test]
        public void FLASHDeconvParameters_DefaultConstructor_AverageResidueModelIsAveragine()
        {
            var p = new FLASHDeconvolutionParameters();
            Assert.That(p.AverageResidueModel, Is.Not.Null);
            Assert.That(p.AverageResidueModel, Is.InstanceOf<Averagine>());
        }

        [Test]
        public void FLASHDeconvParameters_CustomConstructor_SetsAllValues()
        {
            var p = new FLASHDeconvolutionParameters(
                minCharge: 2,
                maxCharge: 30,
                deconvolutionTolerancePpm: 5.0,
                minIsotopicPeakCount: 5,
                maxIsotopicPeakCount: 25,
                precursorHarmonicCount: 3,
                minCosineScore: 0.75,
                isotopeIntensityRatioThreshold: 0.15,
                minMassRange: 200.0,
                maxMassRange: 50000.0,
                polarity: Polarity.Negative);

            Assert.That(p.MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(p.MaxAssumedChargeState, Is.EqualTo(30));
            Assert.That(p.DeconvolutionTolerancePpm, Is.EqualTo(5.0).Within(1e-10));
            Assert.That(p.MinIsotopicPeakCount, Is.EqualTo(5));
            Assert.That(p.MaxIsotopicPeakCount, Is.EqualTo(25));
            Assert.That(p.PrecursorHarmonicCount, Is.EqualTo(3));
            Assert.That(p.MinCosineScore, Is.EqualTo(0.75).Within(1e-10));
            Assert.That(p.IsotopeIntensityRatioThreshold, Is.EqualTo(0.15).Within(1e-10));
            Assert.That(p.MinMassRange, Is.EqualTo(200.0).Within(1e-10));
            Assert.That(p.MaxMassRange, Is.EqualTo(50000.0).Within(1e-10));
            Assert.That(p.Polarity, Is.EqualTo(Polarity.Negative));
        }

        [Test]
        public void FLASHDeconvParameters_InheritsFromDeconvolutionParameters()
        {
            var p = new FLASHDeconvolutionParameters();
            Assert.That(p, Is.InstanceOf<DeconvolutionParameters>());
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. Factory dispatch — Deconvoluter routes to FLASHDeconvAlgorithm
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Deconvoluter_WithFLASHParams_DoesNotThrow()
        {
            var spectrum = SmallPositiveSpectrum();
            var p = new FLASHDeconvolutionParameters();

            IEnumerable<IsotopicEnvelope> result = null;
            Assert.DoesNotThrow(() =>
                result = Deconvoluter.Deconvolute(spectrum, p).ToList());

            Assert.That(result, Is.Not.Null);
        }

        [Test]
        public void Deconvoluter_WithFLASHParams_ReturnsIEnumerableOfIsotopicEnvelope()
        {
            var spectrum = SmallPositiveSpectrum();
            var result = Deconvoluter.Deconvolute(spectrum, new FLASHDeconvolutionParameters());

            Assert.That(result, Is.InstanceOf<IEnumerable<IsotopicEnvelope>>());
        }

        [Test]
        public void Deconvoluter_MsDataScanOverload_WithFLASHParams_DoesNotThrow()
        {
            var spectrum = SmallPositiveSpectrum();
            var scan = new MsDataScan(
                spectrum, 1, 1, true, Polarity.Positive,
                1.0, new MzRange(400, 800), null,
                MZAnalyzerType.Orbitrap, spectrum.SumOfAllY,
                null, null, "scan=1");

            IEnumerable<IsotopicEnvelope> result = null;
            Assert.DoesNotThrow(() =>
                result = Deconvoluter.Deconvolute(scan, new FLASHDeconvolutionParameters()).ToList());

            Assert.That(result, Is.Not.Null);
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. Guard clauses — boundary and edge cases that must never throw
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FLASHDeconv_EmptySpectrum_ReturnsEmptyNotNull()
        {
            var result = Deconvoluter.Deconvolute(
                EmptySpectrum(),
                new FLASHDeconvolutionParameters()).ToList();

            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(0));
        }

        [Test]
        public void FLASHDeconv_NullRange_DefaultsToSpectrumRangeWithoutThrowing()
        {
            var spectrum = SmallPositiveSpectrum();

            IEnumerable<IsotopicEnvelope> result = null;
            Assert.DoesNotThrow(() =>
                result = Deconvoluter.Deconvolute(spectrum, new FLASHDeconvolutionParameters(), null).ToList());

            Assert.That(result, Is.Not.Null);
        }

        [Test]
        public void FLASHDeconv_NegativePolarity_DoesNotThrow()
        {
            var spectrum = SmallPositiveSpectrum();
            var p = new FLASHDeconvolutionParameters(polarity: Polarity.Negative);

            Assert.DoesNotThrow(() =>
                _ = Deconvoluter.Deconvolute(spectrum, p).ToList());
        }

        [Test]
        public void FLASHDeconv_SinglePeakSpectrum_DoesNotThrow()
        {
            var spectrum = new MzSpectrum(
                new double[] { 500.0 },
                new double[] { 1000.0 },
                shouldCopy: false);

            Assert.DoesNotThrow(() =>
                _ = Deconvoluter.Deconvolute(spectrum, new FLASHDeconvolutionParameters()).ToList());
        }

        [Test]
        public void FLASHDeconv_RangeNarrowerThanSpectrum_DoesNotThrow()
        {
            var spectrum = SmallPositiveSpectrum();
            var narrowRange = new MzRange(499.0, 501.5);

            Assert.DoesNotThrow(() =>
                _ = Deconvoluter.Deconvolute(spectrum, new FLASHDeconvolutionParameters(), narrowRange).ToList());
        }

        // ══════════════════════════════════════════════════════════════════════
        // 4. Skeleton behaviour — documents what the skeleton returns today.
        //    Replace assertions in this region when the algorithm is implemented.
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// The skeleton returns empty. This test documents that behaviour and
        /// will need to be updated once the real algorithm is implemented.
        /// </summary>
        [Test]
        public void FLASHDeconv_Skeleton_ReturnsEmptyOnRealSpectrum()
        {
            var result = Deconvoluter.Deconvolute(
                SmallPositiveSpectrum(),
                new FLASHDeconvolutionParameters()).ToList();

            // Skeleton behaviour: always empty.
            // TODO: Replace Is.Empty with real quality assertions after implementation.
            Assert.That(result, Is.Empty,
                "FLASHDeconv skeleton should return empty. " +
                "Update this assertion when the algorithm is implemented.");
        }

        // ══════════════════════════════════════════════════════════════════════
        // 5. Enum — completeness check
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void DeconvolutionType_FLASHDeconvolution_ExistsInEnum()
        {
            var names = System.Enum.GetNames(typeof(DeconvolutionType));
            Assert.That(names, Contains.Item("FLASHDeconvolution"));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 6. NeutralMassSpectrum short-circuit — FLASHDeconv params do not
        //    interfere with the existing neutral-mass pass-through path
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void FLASHDeconv_NeutralMassSpectrum_ShortCircuitsToNeutralMassPath()
        {
            // Arrange: build a NeutralMassSpectrum (pre-deconvoluted)
            var masses = new double[] { 5000.0, 10000.0 };
            var intensities = new double[] { 1000.0, 2000.0 };
            var charges = new int[] { 5, 10 };
            var neutral = new NeutralMassSpectrum(masses, intensities, charges, shouldCopy: false);

            // Act: Deconvoluter should short-circuit to DeconvoluteNeutralMassSpectrum,
            // returning one IsotopicEnvelope per peak regardless of which parameters are used.
            var result = Deconvoluter.Deconvolute(neutral, new FLASHDeconvolutionParameters()).ToList();

            // Assert: both peaks are returned as-is (neutral mass path)
            Assert.That(result.Count, Is.EqualTo(2));
            Assert.That(result[0].MonoisotopicMass, Is.EqualTo(5000.0).Within(1e-6));
            Assert.That(result[1].MonoisotopicMass, Is.EqualTo(10000.0).Within(1e-6));
            Assert.That(result[0].Charge, Is.EqualTo(5));
            Assert.That(result[1].Charge, Is.EqualTo(10));
        }
    }
}