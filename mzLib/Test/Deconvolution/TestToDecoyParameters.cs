using System;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Tests for <see cref="DeconvolutionParameters.ToDecoyParameters"/> and its
    /// implementation in each concrete parameter class.
    ///
    /// The critical behaviours being tested are:
    ///   1. The returned object is a distinct instance of the same concrete type.
    ///   2. All user-configured properties are copied to the decoy instance.
    ///   3. The decoy instance has <see cref="DeconvolutionParameters.ExpectedIsotopeSpacing"/>
    ///      set to <see cref="DecoyAveragine.DefaultDecoyIsotopeSpacing"/> (0.9444 Da).
    ///   4. The original parameters object is not modified.
    ///   5. Repeated calls to ToDecoyParameters() return the same cached instance.
    ///   6. For IsoDecDeconvolutionParameters: the decoy spacing reaches the DLL via
    ///      <see cref="IsoDecDeconvolutionParameters.ToIsoSettings"/>.
    /// </summary>
    [TestFixture]
    public class TestToDecoyParameters
    {
        // ══════════════════════════════════════════════════════════════════════
        // ClassicDeconvolutionParameters
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Classic_ToDecoyParameters_ReturnsDifferentObjectReference()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy, Is.Not.SameAs(original));
        }

        [Test]
        public void Classic_ToDecoyParameters_ReturnsSameConcreteType()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy, Is.InstanceOf<ClassicDeconvolutionParameters>());
        }

        [Test]
        public void Classic_ToDecoyParameters_CopiesAllProperties()
        {
            var original = new ClassicDeconvolutionParameters(
                minCharge: 2,
                maxCharge: 15,
                deconPpm: 7.5,
                intensityRatio: 4.0,
                polarity: Polarity.Negative);

            var decoy = (ClassicDeconvolutionParameters)original.ToDecoyParameters();

            Assert.That(decoy.MinAssumedChargeState, Is.EqualTo(original.MinAssumedChargeState));
            Assert.That(decoy.MaxAssumedChargeState, Is.EqualTo(original.MaxAssumedChargeState));
            Assert.That(decoy.DeconvolutionTolerancePpm, Is.EqualTo(original.DeconvolutionTolerancePpm));
            Assert.That(decoy.IntensityRatioLimit, Is.EqualTo(original.IntensityRatioLimit));
            Assert.That(decoy.Polarity, Is.EqualTo(original.Polarity));
            Assert.That(decoy.DeconvolutionType, Is.EqualTo(original.DeconvolutionType));
        }

        [Test]
        public void Classic_ToDecoyParameters_SetsDecoyIsotopeSpacing()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy.ExpectedIsotopeSpacing,
                Is.EqualTo(DecoyAveragine.DefaultDecoyIsotopeSpacing).Within(1e-10));
        }

        [Test]
        public void Classic_ToDecoyParameters_UsesDecoyAveragineModel()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy.AverageResidueModel, Is.InstanceOf<DecoyAveragine>());
        }

        [Test]
        public void Classic_ToDecoyParameters_DoesNotModifyOriginal()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            _ = original.ToDecoyParameters();

            Assert.That(original.ExpectedIsotopeSpacing,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10),
                "ToDecoyParameters must not modify the original parameters object");
            Assert.That(original.AverageResidueModel, Is.InstanceOf<Averagine>());
        }

        [Test]
        public void Classic_ToDecoyParameters_ReturnsCachedInstance()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy1 = original.ToDecoyParameters();
            var decoy2 = original.ToDecoyParameters();

            Assert.That(decoy1, Is.SameAs(decoy2),
                "Repeated calls to ToDecoyParameters() should return the same cached instance");
        }

        [Test]
        public void Classic_OriginalSpacing_IsC13MinusC12ByDefault()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);

            Assert.That(original.ExpectedIsotopeSpacing,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10));
        }

        // ══════════════════════════════════════════════════════════════════════
        // IsoDecDeconvolutionParameters
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void IsoDec_ToDecoyParameters_ReturnsDifferentObjectReference()
        {
            var original = new IsoDecDeconvolutionParameters();
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy, Is.Not.SameAs(original));
        }

        [Test]
        public void IsoDec_ToDecoyParameters_ReturnsSameConcreteType()
        {
            var original = new IsoDecDeconvolutionParameters();
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy, Is.InstanceOf<IsoDecDeconvolutionParameters>());
        }

        [Test]
        public void IsoDec_ToDecoyParameters_CopiesUserParameters()
        {
            var original = new IsoDecDeconvolutionParameters(
                polarity: Polarity.Positive,
                cssThreshold: 0.85f,
                matchTolerance: 7.0f,
                maxShift: 5,
                knockdownRounds: 3,
                minAreaCovered: 0.30f,
                relativeDataThreshold: 0.10f);

            var decoy = (IsoDecDeconvolutionParameters)original.ToDecoyParameters();

            Assert.That(decoy.CssThreshold, Is.EqualTo(original.CssThreshold));
            Assert.That(decoy.MatchTolerance, Is.EqualTo(original.MatchTolerance));
            Assert.That(decoy.MaxShift, Is.EqualTo(original.MaxShift));
            Assert.That(decoy.KnockdownRounds, Is.EqualTo(original.KnockdownRounds));
            Assert.That(decoy.MinAreaCovered, Is.EqualTo(original.MinAreaCovered));
            Assert.That(decoy.DataThreshold, Is.EqualTo(original.DataThreshold));
            Assert.That(decoy.Polarity, Is.EqualTo(original.Polarity));
            Assert.That(decoy.DeconvolutionType, Is.EqualTo(original.DeconvolutionType));
        }

        [Test]
        public void IsoDec_ToDecoyParameters_SetsDecoyIsotopeSpacing()
        {
            var original = new IsoDecDeconvolutionParameters();
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy.ExpectedIsotopeSpacing,
                Is.EqualTo(DecoyAveragine.DefaultDecoyIsotopeSpacing).Within(1e-10));
        }

        [Test]
        public void IsoDec_ToDecoyParameters_UsesDecoyAveragineModel()
        {
            var original = new IsoDecDeconvolutionParameters();
            var decoy = original.ToDecoyParameters();

            Assert.That(decoy.AverageResidueModel, Is.InstanceOf<DecoyAveragine>());
        }

        [Test]
        public void IsoDec_ToDecoyParameters_DoesNotModifyOriginal()
        {
            var original = new IsoDecDeconvolutionParameters();
            _ = original.ToDecoyParameters();

            Assert.That(original.ExpectedIsotopeSpacing,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10));
            Assert.That(original.AverageResidueModel, Is.InstanceOf<Averagine>());
        }

        [Test]
        public void IsoDec_ToDecoyParameters_ReturnsCachedInstance()
        {
            var original = new IsoDecDeconvolutionParameters();
            var decoy1 = original.ToDecoyParameters();
            var decoy2 = original.ToDecoyParameters();

            Assert.That(decoy1, Is.SameAs(decoy2));
        }

        /// <summary>
        /// The critical test for IsoDecDeconvolutionParameters.
        ///
        /// Verifies that the decoy isotope spacing (0.9444 Da) is correctly forwarded
        /// to the IsoDec DLL via <see cref="IsoDecDeconvolutionParameters.ToIsoSettings"/>.
        /// This works because ToIsoSettings() uses ExpectedIsotopeSpacing rather than
        /// the legacy MassDiffC field for mass_diff_c.
        /// </summary>
        [Test]
        public void IsoDec_ToDecoyParameters_DecoySpacingReachesDll()
        {
            var original = new IsoDecDeconvolutionParameters();

            // Warm the original's cache before cloning. Without this, _isoSettings on
            // the original is null, so any future "cache leaks to decoy" regression
            // (e.g. a MemberwiseClone-based ToDecoyParameters that copies the field)
            // would be invisible because there'd be no cache to copy in the first place.
            var originalSettings = original.ToIsoSettings();
            Assert.That(originalSettings.mass_diff_c, Is.EqualTo(Constants.C13MinusC12).Within(1e-6),
                "Sanity: original cache should hold C13MinusC12 before cloning");

            var decoy = (IsoDecDeconvolutionParameters)original.ToDecoyParameters();
            IsoDecDeconvolutionParameters.IsoSettings decoySettings = decoy.ToIsoSettings();

            Assert.That(decoySettings.mass_diff_c,
                Is.EqualTo(DecoyAveragine.DefaultDecoyIsotopeSpacing).Within(1e-6),
                "The decoy IsoSettings.mass_diff_c must equal 0.9444 Da so the " +
                "shifted spacing reaches the IsoDec DLL");
        }

        [Test]
        public void IsoDec_Original_UsesC13MinusC12InDll()
        {
            var original = new IsoDecDeconvolutionParameters();
            IsoDecDeconvolutionParameters.IsoSettings settings = original.ToIsoSettings();

            Assert.That(settings.mass_diff_c,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-6),
                "Normal (non-decoy) IsoSettings.mass_diff_c should be C13MinusC12");
        }
    }
}