using System;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Unit tests for <see cref="DeconvolutionParameters.ShallowClone"/> and its
    /// override in <see cref="IsoDecDeconvolutionParameters"/>.
    ///
    /// The critical behaviour being tested is that:
    ///   1. The clone is a distinct object of the same concrete type.
    ///   2. All user-configured properties are copied to the clone.
    ///   3. Mutating the clone does not affect the original.
    ///   4. For IsoDecDeconvolutionParameters specifically: the _isoSettings cache
    ///      is nulled on the clone so that ToIsoSettings() rebuilds from the clone's
    ///      fields — in particular from DecoyIsotopeDistance, which MakeDecoyParameters
    ///      sets to 0.9444 after cloning. Without this, the cached struct would still
    ///      hold MassDiffC and the shifted spacing would never reach the DLL.
    /// </summary>
    [TestFixture]
    public class TestShallowClone
    {
        // ══════════════════════════════════════════════════════════════════════
        // ClassicDeconvolutionParameters — exercises the base ShallowClone
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Classic_ShallowClone_ReturnsDifferentObjectReference()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var clone = original.ShallowClone();

            Assert.That(clone, Is.Not.SameAs(original),
                "ShallowClone must return a new object, not the same reference");
        }

        [Test]
        public void Classic_ShallowClone_ReturnsSameConcreteType()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var clone = original.ShallowClone();

            Assert.That(clone, Is.InstanceOf<ClassicDeconvolutionParameters>(),
                "ShallowClone must return the same concrete type, not just the base class");
        }

        [Test]
        public void Classic_ShallowClone_CopiesAllProperties()
        {
            var original = new ClassicDeconvolutionParameters(
                minCharge: 2,
                maxCharge: 15,
                deconPpm: 7.5,
                intensityRatio: 4.0,
                polarity: Polarity.Negative);

            var clone = (ClassicDeconvolutionParameters)original.ShallowClone();

            Assert.That(clone.MinAssumedChargeState, Is.EqualTo(original.MinAssumedChargeState));
            Assert.That(clone.MaxAssumedChargeState, Is.EqualTo(original.MaxAssumedChargeState));
            Assert.That(clone.DeconvolutionTolerancePpm, Is.EqualTo(original.DeconvolutionTolerancePpm));
            Assert.That(clone.IntensityRatioLimit, Is.EqualTo(original.IntensityRatioLimit));
            Assert.That(clone.Polarity, Is.EqualTo(original.Polarity));
            Assert.That(clone.DeconvolutionType, Is.EqualTo(original.DeconvolutionType));
        }

        [Test]
        public void Classic_ShallowClone_DefaultDecoyIsotopeDistance_EqualsC13MinusC12()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var clone = original.ShallowClone();

            Assert.That(clone.DecoyIsotopeDistance,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10),
                "A fresh clone should carry the default DecoyIsotopeDistance (C13MinusC12), " +
                "not the decoy value");
        }

        [Test]
        public void Classic_ShallowClone_MutatingClone_DoesNotAffectOriginal()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var clone = original.ShallowClone();

            // Mutate the clone via the decoy generator path
            clone.DecoyIsotopeDistance = 0.9444;
            clone.IsDecoyRun = true;

            // Original must be unchanged
            Assert.That(original.DecoyIsotopeDistance,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10),
                "Mutating the clone's DecoyIsotopeDistance must not affect the original");
            Assert.That(original.IsDecoyRun, Is.False,
                "Mutating the clone's IsDecoyRun must not affect the original");
        }

        // ══════════════════════════════════════════════════════════════════════
        // IsoDecDeconvolutionParameters — exercises the override
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void IsoDec_ShallowClone_ReturnsDifferentObjectReference()
        {
            var original = new IsoDecDeconvolutionParameters();
            var clone = original.ShallowClone();

            Assert.That(clone, Is.Not.SameAs(original));
        }

        [Test]
        public void IsoDec_ShallowClone_ReturnsSameConcreteType()
        {
            var original = new IsoDecDeconvolutionParameters();
            var clone = original.ShallowClone();

            Assert.That(clone, Is.InstanceOf<IsoDecDeconvolutionParameters>(),
                "ShallowClone must return IsoDecDeconvolutionParameters, not the base type");
        }

        [Test]
        public void IsoDec_ShallowClone_CopiesUserParameters()
        {
            var original = new IsoDecDeconvolutionParameters(
                polarity: Polarity.Positive,
                cssThreshold: 0.85f,
                matchTolerance: 7.0f,
                maxShift: 5,
                knockdownRounds: 3,
                minAreaCovered: 0.30f,
                relativeDataThreshold: 0.10f);

            var clone = (IsoDecDeconvolutionParameters)original.ShallowClone();

            Assert.That(clone.CssThreshold, Is.EqualTo(original.CssThreshold));
            Assert.That(clone.MatchTolerance, Is.EqualTo(original.MatchTolerance));
            Assert.That(clone.MaxShift, Is.EqualTo(original.MaxShift));
            Assert.That(clone.KnockdownRounds, Is.EqualTo(original.KnockdownRounds));
            Assert.That(clone.MinAreaCovered, Is.EqualTo(original.MinAreaCovered));
            Assert.That(clone.DataThreshold, Is.EqualTo(original.DataThreshold));
            Assert.That(clone.Polarity, Is.EqualTo(original.Polarity));
            Assert.That(clone.DeconvolutionType, Is.EqualTo(original.DeconvolutionType));
        }

        [Test]
        public void IsoDec_ShallowClone_CopiesHardCodedParameters()
        {
            var original = new IsoDecDeconvolutionParameters();
            var clone = (IsoDecDeconvolutionParameters)original.ShallowClone();

            // Hard-coded parameters should also be present on the clone
            Assert.That(clone.MassDiffC, Is.EqualTo(original.MassDiffC));
            Assert.That(clone.IsoLength, Is.EqualTo(original.IsoLength));
            Assert.That(clone.AdductMass, Is.EqualTo(original.AdductMass).Within(1e-6f));
            Assert.That(clone.MinPeaks, Is.EqualTo(original.MinPeaks));
        }

        /// <summary>
        /// This is the critical test for the IsoDecDeconvolutionParameters override.
        ///
        /// The purpose of overriding ShallowClone is to null the _isoSettings cache
        /// on the clone. If the cache is NOT nulled, ToIsoSettings() on the clone will
        /// return the cached struct built from the original's MassDiffC, ignoring any
        /// subsequent change to DecoyIsotopeDistance. This would mean the decoy isotope
        /// spacing never reaches the IsoDec DLL.
        ///
        /// We verify this by:
        ///   1. Priming the cache on the original by calling ToIsoSettings().
        ///   2. Cloning the original.
        ///   3. Setting DecoyIsotopeDistance = 0.9444 on the clone.
        ///   4. Calling ToIsoSettings() on the clone.
        ///   5. Asserting that the clone's IsoSettings.mass_diff_c == 0.9444,
        ///      NOT the original MassDiffC value.
        /// </summary>
        [Test]
        public void IsoDec_ShallowClone_NullsIsoSettingsCache_SoDecoyDistanceReachesDll()
        {
            var original = new IsoDecDeconvolutionParameters();

            // Prime the cache on the original.
            // ToIsoSettings() uses DecoyIsotopeDistance (which defaults to
            // Constants.C13MinusC12 = 1.003354838...), NOT the rounded MassDiffC
            // property (1.0033). These are intentionally different values.
            IsoDecDeconvolutionParameters.IsoSettings originalSettings = original.ToIsoSettings();
            Assert.That(originalSettings.mass_diff_c,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-6),
                "Precondition: original ToIsoSettings() should use DecoyIsotopeDistance " +
                "(= Constants.C13MinusC12 by default), not the rounded MassDiffC property");

            // Clone and apply the decoy isotope distance (mimicking MakeDecoyParameters)
            const double decoyDistance = 0.9444;
            var clone = (IsoDecDeconvolutionParameters)original.ShallowClone();
            clone.DecoyIsotopeDistance = decoyDistance;

            // ToIsoSettings() on the clone must rebuild from DecoyIsotopeDistance
            IsoDecDeconvolutionParameters.IsoSettings cloneSettings = clone.ToIsoSettings();

            Assert.That(cloneSettings.mass_diff_c,
                Is.EqualTo(decoyDistance).Within(1e-6),
                "After ShallowClone + DecoyIsotopeDistance = 0.9444, ToIsoSettings() " +
                "on the clone must return mass_diff_c = 0.9444. If it returns " +
                $"{original.MassDiffC}, the _isoSettings cache was not nulled.");
        }

        [Test]
        public void IsoDec_ShallowClone_OriginalCacheUnaffectedAfterCloningAndMutating()
        {
            var original = new IsoDecDeconvolutionParameters();

            // Prime the original's cache
            IsoDecDeconvolutionParameters.IsoSettings beforeClone = original.ToIsoSettings();

            // Clone and mutate the clone's isotope distance
            var clone = (IsoDecDeconvolutionParameters)original.ShallowClone();
            clone.DecoyIsotopeDistance = 0.9444;
            _ = clone.ToIsoSettings(); // force clone to rebuild its settings

            // Original's ToIsoSettings() must still return the original MassDiffC
            IsoDecDeconvolutionParameters.IsoSettings afterClone = original.ToIsoSettings();

            Assert.That(afterClone.mass_diff_c,
                Is.EqualTo(beforeClone.mass_diff_c).Within(1e-10),
                "Mutating the clone and rebuilding its IsoSettings must not affect " +
                "the original's cached IsoSettings");
        }

        [Test]
        public void IsoDec_ShallowClone_MutatingClone_DoesNotAffectOriginalBaseProperties()
        {
            var original = new IsoDecDeconvolutionParameters();
            var clone = original.ShallowClone();

            clone.DecoyIsotopeDistance = 0.9444;
            clone.IsDecoyRun = true;

            Assert.That(original.DecoyIsotopeDistance,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10));
            Assert.That(original.IsDecoyRun, Is.False);
        }

        // ══════════════════════════════════════════════════════════════════════
        // MakeDecoyParameters — exercises the full path through ShallowClone
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void MakeDecoyParameters_Classic_ReturnsSameConcreteType()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy = DeconvolutionDecoyGenerator.MakeDecoyParameters(original);

            Assert.That(decoy, Is.InstanceOf<ClassicDeconvolutionParameters>());
        }

        [Test]
        public void MakeDecoyParameters_Classic_SetsDecoyIsotopeDistanceAndFlag()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var decoy = DeconvolutionDecoyGenerator.MakeDecoyParameters(original);

            Assert.That(decoy.DecoyIsotopeDistance,
                Is.EqualTo(DeconvolutionDecoyGenerator.DecoyIsotopeDistance).Within(1e-10));
            Assert.That(decoy.IsDecoyRun, Is.True);
        }

        [Test]
        public void MakeDecoyParameters_Classic_DoesNotModifyOriginal()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            _ = DeconvolutionDecoyGenerator.MakeDecoyParameters(original);

            Assert.That(original.DecoyIsotopeDistance,
                Is.EqualTo(Constants.C13MinusC12).Within(1e-10),
                "MakeDecoyParameters must not modify the original parameters object");
            Assert.That(original.IsDecoyRun, Is.False);
        }

        [Test]
        public void MakeDecoyParameters_IsoDec_ReturnsSameConcreteType()
        {
            var original = new IsoDecDeconvolutionParameters();
            var decoy = DeconvolutionDecoyGenerator.MakeDecoyParameters(original);

            Assert.That(decoy, Is.InstanceOf<IsoDecDeconvolutionParameters>());
        }

        [Test]
        public void MakeDecoyParameters_IsoDec_DecoyDistanceReachesDll()
        {
            var original = new IsoDecDeconvolutionParameters();

            // Prime the original's cache before cloning
            _ = original.ToIsoSettings();

            var decoy = (IsoDecDeconvolutionParameters)
                DeconvolutionDecoyGenerator.MakeDecoyParameters(original);

            IsoDecDeconvolutionParameters.IsoSettings decoySettings = decoy.ToIsoSettings();

            Assert.That(decoySettings.mass_diff_c,
                Is.EqualTo(DeconvolutionDecoyGenerator.DecoyIsotopeDistance).Within(1e-6),
                "After MakeDecoyParameters, the decoy IsoSettings.mass_diff_c must " +
                "equal 0.9444, not the original MassDiffC");
        }
    }
}