using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Tests for <see cref="MultipleDeconParameters"/> and <see cref="MultipleDeconvolutionAlgorithm"/>.
    ///
    /// Coverage:
    ///   1. <see cref="MultipleDeconParameters"/> construction and property correctness.
    ///   2. <see cref="DeconvolutionType.Multiple"/> is dispatched to <see cref="MultipleDeconvolutionAlgorithm"/>.
    ///   3. A single algorithm wrapped in Multiple produces the same envelopes as calling it directly.
    ///   4. Results from all wrapped algorithms are concatenated (Classic + IsoDec).
    ///   5. The MzRange filter is forwarded correctly to each sub-algorithm.
    ///   6. The <see cref="MsDataScan"/> overload works end-to-end.
    ///   7. A FromFile sub-parameter inside Multiple triggers RT-range synthesis in
    ///      <see cref="Deconvoluter.Deconvolute(MsDataScan, DeconvolutionParameters, MzRange)"/>.
    ///   8. <see cref="MultipleDeconParameters.ToDecoyParameters"/> returns a non-null
    ///      <see cref="MultipleDeconParameters"/> when every sub-param supports decoys.
    ///   9. <see cref="MultipleDeconParameters.ToDecoyParameters"/> returns <c>null</c>
    ///      when any sub-param does NOT support decoys (e.g. FromFile).
    ///  10. The decoy <see cref="MultipleDeconParameters"/> carries the correct
    ///      <see cref="DecoyAveragine.DefaultDecoyIsotopeSpacing"/> on every sub-param.
    ///  11. Charge and polarity are copied faithfully to the decoy instance.
    ///  12. Repeated calls to <c>ToDecoyParameters()</c> return the same cached instance.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestMultipleDeconvolution
    {
        // ?? shared helpers ????????????????????????????????????????????????????

        private static MzSpectrum LoadProteoformSpectrum()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"DataFiles\14kDaProteoformMzIntensityMs1.txt");
            string[] lines = File.ReadAllLines(path);
            double[] mzs = new double[lines.Length];
            double[] intensities = new double[lines.Length];
            for (int i = 0; i < lines.Length; i++)
            {
                string[] parts = lines[i].Split('\t');
                mzs[i] = double.Parse(parts[0], CultureInfo.InvariantCulture);
                intensities[i] = double.Parse(parts[1], CultureInfo.InvariantCulture);
            }
            return new MzSpectrum(mzs, intensities, false);
        }

        private static MsDataScan MakeProteoformScan(MzSpectrum spectrum)
            => new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0,
                new MzRange(495, 1617), "scan", MZAnalyzerType.Unknown,
                spectrum.SumOfAllY, null, null, null,
                740.372202090153, 19, 108419280, 740.37, 4);

        private static ClassicDeconvolutionParameters DefaultClassicParams()
            => new ClassicDeconvolutionParameters(1, 60, 4, 3);

        private static IsoDecDeconvolutionParameters DefaultIsoDecParams()
            => new IsoDecDeconvolutionParameters();

        // ?? 1. Construction ???????????????????????????????????????????????????

        [Test]
        public void MultipleDeconParameters_DeconvolutionType_IsMultiple()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            Assert.That(p.DeconvolutionType, Is.EqualTo(DeconvolutionType.Multiple));
        }

        [Test]
        public void MultipleDeconParameters_WrapsAllProvidedParameters()
        {
            var classic = DefaultClassicParams();
            var isodec = DefaultIsoDecParams();

            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classic, isodec }, 1, 60);

            Assert.That(p.Parameters.Length, Is.EqualTo(2));
            Assert.That(p.Parameters[0], Is.SameAs(classic));
            Assert.That(p.Parameters[1], Is.SameAs(isodec));
        }

        [Test]
        public void MultipleDeconParameters_ChargeAndPolarity_StoredCorrectly()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() },
                minCharge: 2, maxCharge: 30, polarity: Polarity.Negative);

            Assert.That(p.MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(p.MaxAssumedChargeState, Is.EqualTo(30));
            Assert.That(p.Polarity, Is.EqualTo(Polarity.Negative));
        }

        // ?? 2. Algorithm dispatch ?????????????????????????????????????????????

        [Test]
        public void MultipleDecon_Dispatch_DoesNotThrow()
        {
            var spectrum = LoadProteoformSpectrum();
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            Assert.DoesNotThrow(() => _ = Deconvoluter.Deconvolute(spectrum, p).ToList());
        }

        // ?? 3. Single wrapped algorithm matches direct call ???????????????????

        [Test]
        public void MultipleDecon_SingleWrappedClassic_MatchesDirectClassic()
        {
            var spectrum = LoadProteoformSpectrum();
            var classicParams = DefaultClassicParams();
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classicParams }, 1, 60);

            var direct = Deconvoluter.Deconvolute(spectrum, classicParams).ToList();
            var wrapped = Deconvoluter.Deconvolute(spectrum, multipleParams).ToList();

            // Every envelope produced by the direct call should appear in the wrapped result.
            Assert.That(wrapped.Count, Is.EqualTo(direct.Count));
            Assert.That(
                direct.Select(e => e.MonoisotopicMass).OrderBy(m => m),
                Is.EquivalentTo(wrapped.Select(e => e.MonoisotopicMass).OrderBy(m => m)));
        }

        [Test]
        public void MultipleDecon_SingleWrappedIsoDec_MatchesDirectIsoDec()
        {
            var spectrum = LoadProteoformSpectrum();
            var isoDecParams = DefaultIsoDecParams();
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { isoDecParams }, 1, 60);

            var direct = Deconvoluter.Deconvolute(spectrum, isoDecParams).ToList();
            var wrapped = Deconvoluter.Deconvolute(spectrum, multipleParams).ToList();

            Assert.That(wrapped.Count, Is.EqualTo(direct.Count));
        }

        // ?? 4. Results from all sub-algorithms are concatenated ???????????????

        [Test]
        public void MultipleDecon_ClassicAndIsoDec_CountIsAtLeastSumOfIndividual()
        {
            // The combined result may contain duplicate envelopes (the same peak found
            // by both algorithms) but must contain at least as many as each alone.
            var spectrum = LoadProteoformSpectrum();
            var classicParams = DefaultClassicParams();
            var isoDecParams = DefaultIsoDecParams();
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classicParams, isoDecParams }, 1, 60);

            int classicCount = Deconvoluter.Deconvolute(spectrum, classicParams).Count();
            int isoDecCount = Deconvoluter.Deconvolute(spectrum, isoDecParams).Count();
            int multipleCount = Deconvoluter.Deconvolute(spectrum, multipleParams).Count();

            Assert.That(multipleCount, Is.EqualTo(classicCount + isoDecCount));
        }

        [Test]
        public void MultipleDecon_ThreeAlgorithms_ConcatenatesAll()
        {
            var spectrum = LoadProteoformSpectrum();
            var classic1 = DefaultClassicParams();
            var classic2 = new ClassicDeconvolutionParameters(1, 30, 4, 3);
            var isoDecParams = DefaultIsoDecParams();

            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classic1, classic2, isoDecParams }, 1, 60);

            int expected = Deconvoluter.Deconvolute(spectrum, classic1).Count()
                         + Deconvoluter.Deconvolute(spectrum, classic2).Count()
                         + Deconvoluter.Deconvolute(spectrum, isoDecParams).Count();

            int actual = Deconvoluter.Deconvolute(spectrum, multipleParams).Count();

            Assert.That(actual, Is.EqualTo(expected));
        }

        // ?? 5. MzRange is forwarded correctly ?????????????????????????????????

        [Test]
        public void MultipleDecon_WithMzRange_RespectsRange()
        {
            var spectrum = LoadProteoformSpectrum();
            var narrowRange = new MzRange(700, 800);
            var classicParams = DefaultClassicParams();
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classicParams }, 1, 60);

            var directNarrow = Deconvoluter.Deconvolute(spectrum, classicParams, narrowRange).ToList();
            var wrappedNarrow = Deconvoluter.Deconvolute(spectrum, multipleParams, narrowRange).ToList();

            Assert.That(wrappedNarrow.Count, Is.EqualTo(directNarrow.Count));
            // Every returned peak must fall within the narrow range.
            foreach (var env in wrappedNarrow)
                Assert.That(env.Peaks.All(pk => narrowRange.Contains(pk.mz)));
        }

        [Test]
        public void MultipleDecon_FullRangeVsNarrowRange_NarrowHasFewerResults()
        {
            var spectrum = LoadProteoformSpectrum();
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            int fullCount = Deconvoluter.Deconvolute(spectrum, multipleParams).Count();
            int narrowCount = Deconvoluter.Deconvolute(spectrum, multipleParams,
                new MzRange(700, 800)).Count();

            Assert.That(narrowCount, Is.LessThan(fullCount));
        }

        // ?? 6. MsDataScan overload ?????????????????????????????????????????????

        [Test]
        public void MultipleDecon_ScanOverload_MatchesSpectrumOverload()
        {
            var spectrum = LoadProteoformSpectrum();
            var scan = MakeProteoformScan(spectrum);
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            var fromScan = Deconvoluter.Deconvolute(scan, multipleParams).ToList();
            var fromSpectrum = Deconvoluter.Deconvolute(spectrum, multipleParams).ToList();

            Assert.That(fromScan.Count, Is.EqualTo(fromSpectrum.Count));
        }

        [Test]
        public void MultipleDecon_ScanOverload_NegativeMode_ReturnsResults()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DataFiles", "GUACUG_NegativeMode_Sliced.mzML");
            var scan = MsDataFileReader.GetDataFile(filePath).GetAllScansList().First();

            var classicNeg = new ClassicDeconvolutionParameters(-10, -1, 20, 3, Polarity.Negative);
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classicNeg }, -10, -1, Polarity.Negative);

            var results = Deconvoluter.Deconvolute(scan, multipleParams).ToList();

            Assert.That(results, Is.Not.Empty);
        }

        // ?? 7. FromFile sub-param triggers RT-range synthesis ?????????????????

        [Test]
        public void MultipleDecon_ContainingFromFile_DoesNotThrow_WhenCalledWithScan()
        {
            // Build a minimal synthetic FromFileDeconvolutionParameters (no real file needed).
            var fromFileParams = new FromFileDeconvolutionParameters(
                Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var classicParams = DefaultClassicParams();
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { classicParams, fromFileParams }, 1, 60);

            var spectrum = LoadProteoformSpectrum();
            var scan = MakeProteoformScan(spectrum);

            // The scan overload must synthesize an MzRtRange and NOT throw.
            Assert.DoesNotThrow(() => _ = Deconvoluter.Deconvolute(scan, multipleParams).ToList());
        }

        [Test]
        public void MultipleDecon_ContainingFromFile_WithoutScan_ThrowsArgumentException()
        {
            // Calling the MzSpectrum overload with a FromFile sub-param (and no MzRtRange)
            // must surface the ArgumentException from Deconvoluter.
            var fromFileParams = new FromFileDeconvolutionParameters(
                Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var multipleParams = new MultipleDeconParameters(
                new DeconvolutionParameters[] { fromFileParams }, 1, 60);

            var spectrum = LoadProteoformSpectrum();

            // MultipleDeconvolutionAlgorithm delegates each sub-param back to
            // Deconvoluter.Deconvolute(spectrum, ...), which throws for FromFile
            // without an MzRtRange.
            Assert.Throws<ArgumentException>(
                () => _ = Deconvoluter.Deconvolute(spectrum, multipleParams).ToList());
        }

        // ?? 8. ToDecoyParameters – all sub-params support decoys ??????????????

        [Test]
        public void MultipleDecon_ToDecoyParameters_AllSupport_ReturnsNonNull()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams(), DefaultIsoDecParams() }, 1, 60);

            var decoy = p.ToDecoyParameters();

            Assert.That(decoy, Is.Not.Null);
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_AllSupport_ReturnsMultipleDeconParameters()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams(), DefaultIsoDecParams() }, 1, 60);

            var decoy = p.ToDecoyParameters();

            Assert.That(decoy, Is.InstanceOf<MultipleDeconParameters>());
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_AllSupport_ContainsSameNumberOfSubParams()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams(), DefaultIsoDecParams() }, 1, 60);

            var decoy = (MultipleDeconParameters)p.ToDecoyParameters()!;

            Assert.That(decoy.Parameters.Length, Is.EqualTo(p.Parameters.Length));
        }

        // ?? 9. ToDecoyParameters – any sub-param lacks decoy support ??????????

        [Test]
        public void MultipleDecon_ToDecoyParameters_FromFileSubParam_ReturnsNull()
        {
            // FromFileDeconvolutionParameters.ToDecoyParameters() always returns null.
            var fromFile = new FromFileDeconvolutionParameters(
                Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams(), fromFile }, 1, 60);

            var decoy = p.ToDecoyParameters();

            Assert.That(decoy, Is.Null);
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_AllFromFile_ReturnsNull()
        {
            var ff1 = new FromFileDeconvolutionParameters(Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var ff2 = new FromFileDeconvolutionParameters(Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var p = new MultipleDeconParameters(new DeconvolutionParameters[] { ff1, ff2 }, 1, 60);

            Assert.That(p.ToDecoyParameters(), Is.Null);
        }

        // ?? 10. Decoy isotope spacing ?????????????????????????????????????????

        [Test]
        public void MultipleDecon_ToDecoyParameters_SubParamsHaveDecoyIsotopeSpacing()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams(), DefaultIsoDecParams() }, 1, 60);

            var decoy = (MultipleDeconParameters)p.ToDecoyParameters()!;

            foreach (var sub in decoy.Parameters)
                Assert.That(sub.ExpectedIsotopeSpacing,
                    Is.EqualTo(DecoyAveragine.DefaultDecoyIsotopeSpacing).Within(1e-10));
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_SubParamsUseDecoyAveragineModel()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            var decoy = (MultipleDeconParameters)p.ToDecoyParameters()!;

            foreach (var sub in decoy.Parameters)
                Assert.That(sub.AverageResidueModel, Is.InstanceOf<DecoyAveragine>());
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_TopLevelHasDecoyIsotopeSpacing()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            var decoy = p.ToDecoyParameters()!;

            Assert.That(decoy.ExpectedIsotopeSpacing,
                Is.EqualTo(DecoyAveragine.DefaultDecoyIsotopeSpacing).Within(1e-10));
        }

        // ?? 11. Charge and polarity copied faithfully ?????????????????????????

        [Test]
        public void MultipleDecon_ToDecoyParameters_CopiesChargeStateAndPolarity()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { new ClassicDeconvolutionParameters(-10, -1, 20, 3, Polarity.Negative) },
                minCharge: -10, maxCharge: -1, polarity: Polarity.Negative);

            var decoy = p.ToDecoyParameters()!;

            Assert.That(decoy.MinAssumedChargeState, Is.EqualTo(p.MinAssumedChargeState));
            Assert.That(decoy.MaxAssumedChargeState, Is.EqualTo(p.MaxAssumedChargeState));
            Assert.That(decoy.Polarity, Is.EqualTo(p.Polarity));
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_DoesNotModifyOriginal()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            _ = p.ToDecoyParameters();

            Assert.That(p.ExpectedIsotopeSpacing,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10));
            foreach (var sub in p.Parameters)
                Assert.That(sub.ExpectedIsotopeSpacing,
                    Is.EqualTo(Constants.C13MinusC12).Within(1e-10));
        }

        // ?? 12. Cached instance ???????????????????????????????????????????????

        [Test]
        public void MultipleDecon_ToDecoyParameters_ReturnsCachedInstance()
        {
            var p = new MultipleDeconParameters(
                new DeconvolutionParameters[] { DefaultClassicParams() }, 1, 60);

            var decoy1 = p.ToDecoyParameters();
            var decoy2 = p.ToDecoyParameters();

            Assert.That(decoy1, Is.SameAs(decoy2));
        }
    }
}
