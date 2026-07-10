using System;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test
{
    [TestFixture]
    public class TestDeconvolutionParametersEqualityAndClone
    {
        // ══════════════════════════════════════════════════════════════════════
        // ClassicDeconvolutionParameters
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void Classic_Equal_SameValues_AreEqual()
        {
            var a = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var b = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            Assert.That(a, Is.EqualTo(b));
        }

        [Test]
        public void Classic_Equal_DifferentIntensityRatio_AreNotEqual()
        {
            var a = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var b = new ClassicDeconvolutionParameters(1, 12, 4.0, 5.0);
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void Classic_Equal_DifferentTolerancePpm_AreNotEqual()
        {
            var a = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var b = new ClassicDeconvolutionParameters(1, 12, 7.0, 3.0);
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void Classic_Equal_SameReference_IsEqual()
        {
            var a = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            Assert.That(a, Is.EqualTo(a));
        }

        [Test]
        public void Classic_Equal_Null_IsNotEqual()
        {
            var a = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            Assert.That(a.Equals(null), Is.False);
        }

        [Test]
        public void Classic_Equal_CrossType_IsNotEqual()
        {
            var classic = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var isoDec = new IsoDecDeconvolutionParameters();
            Assert.That(classic, Is.Not.EqualTo(isoDec));
            Assert.That(isoDec, Is.Not.EqualTo(classic));
        }

        [Test]
        public void Classic_GetHashCode_EqualObjects_Match()
        {
            var a = new ClassicDeconvolutionParameters(2, 15, 7.5, 4.0, Polarity.Negative);
            var b = new ClassicDeconvolutionParameters(2, 15, 7.5, 4.0, Polarity.Negative);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void Classic_Clone_ReturnsNewReference()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var clone = original.Clone();
            Assert.That(clone, Is.Not.SameAs(original));
        }

        [Test]
        public void Classic_Clone_EqualsOriginal()
        {
            var original = new ClassicDeconvolutionParameters(2, 15, 7.5, 4.0, Polarity.Negative);
            var clone = original.Clone();
            Assert.That(clone, Is.EqualTo(original));
        }

        [Test]
        public void Classic_Clone_ModificationDoesNotAffectOriginal()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            var clone = original.Clone();
            clone.MinAssumedChargeState = 99;
            Assert.That(original.MinAssumedChargeState, Is.EqualTo(1));
        }

        [Test]
        public void Classic_Clone_PreservesUseGenericScore()
        {
            var original = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0)
            {
                UseGenericScore = true
            };
            var clone = original.Clone();
            Assert.That(clone.UseGenericScore, Is.True);
        }

        // ══════════════════════════════════════════════════════════════════════
        // IsoDecDeconvolutionParameters
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void IsoDec_Equal_SameValues_AreEqual()
        {
            var a = new IsoDecDeconvolutionParameters(
                polarity: Polarity.Positive,
                cssThreshold: 0.85f,
                matchTolerance: 7.0f,
                maxShift: 5);
            var b = new IsoDecDeconvolutionParameters(
                polarity: Polarity.Positive,
                cssThreshold: 0.85f,
                matchTolerance: 7.0f,
                maxShift: 5);
            Assert.That(a, Is.EqualTo(b));
        }

        [Test]
        public void IsoDec_Equal_DifferentCssThreshold_AreNotEqual()
        {
            var a = new IsoDecDeconvolutionParameters(polarity: Polarity.Positive, cssThreshold: 0.7f);
            var b = new IsoDecDeconvolutionParameters(polarity: Polarity.Positive, cssThreshold: 0.85f);
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void IsoDec_Equal_SameReference_IsEqual()
        {
            var a = new IsoDecDeconvolutionParameters();
            Assert.That(a, Is.EqualTo(a));
        }

        [Test]
        public void IsoDec_Equal_Null_IsNotEqual()
        {
            var a = new IsoDecDeconvolutionParameters();
            Assert.That(a.Equals(null), Is.False);
        }

        [Test]
        public void IsoDec_GetHashCode_EqualObjects_Match()
        {
            var a = new IsoDecDeconvolutionParameters(polarity: Polarity.Positive, cssThreshold: 0.8f);
            var b = new IsoDecDeconvolutionParameters(polarity: Polarity.Positive, cssThreshold: 0.8f);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void IsoDec_Clone_ReturnsNewReference()
        {
            var original = new IsoDecDeconvolutionParameters();
            var clone = original.Clone();
            Assert.That(clone, Is.Not.SameAs(original));
        }

        [Test]
        public void IsoDec_Clone_EqualsOriginal()
        {
            var original = new IsoDecDeconvolutionParameters(
                polarity: Polarity.Positive,
                cssThreshold: 0.85f,
                matchTolerance: 7.0f,
                maxShift: 5,
                knockdownRounds: 3,
                minAreaCovered: 0.30f,
                relativeDataThreshold: 0.10f);
            var clone = original.Clone();
            Assert.That(clone, Is.EqualTo(original));
        }

        [Test]
        public void IsoDec_Clone_MzWindowIsDeepCopy()
        {
            var original = new IsoDecDeconvolutionParameters();
            var clone = original.Clone();
            clone.MzWindow[0] = 9999f;
            Assert.That(original.MzWindow[0], Is.Not.EqualTo(9999f));
        }

        [Test]
        public void IsoDec_Clone_PreservesUseGenericScore()
        {
            var original = new IsoDecDeconvolutionParameters
            {
                UseGenericScore = true
            };
            var clone = original.Clone();
            Assert.That(clone.UseGenericScore, Is.True);
        }

        // ══════════════════════════════════════════════════════════════════════
        // MultipleDeconParameters
        // ══════════════════════════════════════════════════════════════════════

        private static ClassicDeconvolutionParameters DefaultClassic()
            => new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);

        private static ClassicDeconvolutionParameters OtherClassic()
            => new ClassicDeconvolutionParameters(1, 12, 7.0, 5.0);

        [Test]
        public void MultipleDecon_Equal_SameValues_AreEqual()
        {
            var sub = DefaultClassic();
            var a = new MultipleDeconParameters(new[] { sub }, 1, 12);
            var b = new MultipleDeconParameters(new[] { sub.Clone() }, 1, 12);
            Assert.That(a, Is.EqualTo(b));
        }

        [Test]
        public void MultipleDecon_Equal_DifferentSubParams_AreNotEqual()
        {
            var a = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12);
            var b = new MultipleDeconParameters(new[] { OtherClassic() }, 1, 12);
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void MultipleDecon_Equal_DifferentSubCount_AreNotEqual()
        {
            var a = new MultipleDeconParameters(new[] { DefaultClassic(), OtherClassic() }, 1, 12);
            var b = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12);
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void MultipleDecon_Equal_SameReference_IsEqual()
        {
            var a = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12);
            Assert.That(a, Is.EqualTo(a));
        }

        [Test]
        public void MultipleDecon_Equal_Null_IsNotEqual()
        {
            var a = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12);
            Assert.That(a.Equals(null), Is.False);
        }

        [Test]
        public void MultipleDecon_Clone_ReturnsNewReference()
        {
            var original = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12);
            var clone = original.Clone();
            Assert.That(clone, Is.Not.SameAs(original));
        }

        [Test]
        public void MultipleDecon_Clone_EqualsOriginal()
        {
            var original = new MultipleDeconParameters(new[] { DefaultClassic(), OtherClassic() }, 1, 12);
            var clone = original.Clone();
            Assert.That(clone, Is.EqualTo(original));
        }

        [Test]
        public void MultipleDecon_Clone_SubParamsAreDeepCopy()
        {
            var original = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12);
            var clone = original.Clone();
            clone.Parameters[0].MinAssumedChargeState = 99;
            Assert.That(original.Parameters[0].MinAssumedChargeState, Is.EqualTo(1));
        }

        [Test]
        public void MultipleDecon_Clone_PreservesUseGenericScore()
        {
            var original = new MultipleDeconParameters(new[] { DefaultClassic() }, 1, 12)
            {
                UseGenericScore = true
            };
            var clone = original.Clone();
            Assert.That(clone.UseGenericScore, Is.True);
        }

        [Test]
        public void MultipleDecon_ToDecoyParameters_AllFromFile_ReturnsNull()
        {
            var ff1 = new FromFileDeconvolutionParameters(Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var ff2 = new FromFileDeconvolutionParameters(Enumerable.Empty<ISingleChargeMs1Feature>(), 1, 60);
            var p = new MultipleDeconParameters(new DeconvolutionParameters[] { ff1, ff2 }, 1, 60);
            Assert.That(p.ToDecoyParameters(), Is.Null);
        }

        // ══════════════════════════════════════════════════════════════════════
        // ExampleNewDeconvolutionParametersTemplate
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void ExampleTemplate_Equal_SameValues_AreEqual()
        {
            var a = new ExampleNewDeconvolutionParametersTemplate(1, 12);
            var b = new ExampleNewDeconvolutionParametersTemplate(1, 12);
            Assert.That(a, Is.EqualTo(b));
        }

        [Test]
        public void ExampleTemplate_Equal_DifferentMinCharge_AreNotEqual()
        {
            var a = new ExampleNewDeconvolutionParametersTemplate(1, 12);
            var b = new ExampleNewDeconvolutionParametersTemplate(2, 12);
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void ExampleTemplate_Clone_ReturnsNewReference()
        {
            var original = new ExampleNewDeconvolutionParametersTemplate(1, 12);
            var clone = original.Clone();
            Assert.That(clone, Is.Not.SameAs(original));
        }

        [Test]
        public void ExampleTemplate_Clone_EqualsOriginal()
        {
            var original = new ExampleNewDeconvolutionParametersTemplate(2, 15, Polarity.Negative)
            {
                UseGenericScore = true
            };
            var clone = original.Clone();
            Assert.That(clone, Is.EqualTo(original));
        }

        // ══════════════════════════════════════════════════════════════════════
        // FromFileDeconvolutionParameters
        // ══════════════════════════════════════════════════════════════════════

        private static SingleChargeMs1Feature Feature(double mz = 600.0, int charge = 2,
            double rtStart = 10.0, double rtEnd = 15.0, double intensity = 1e5)
            => new SingleChargeMs1Feature(mz, charge, rtStart, rtEnd, intensity);

        [Test]
        public void FromFile_Equal_SameFeatures_AreEqual()
        {
            var feats = new[] { Feature(), Feature(601.0, 3) };
            var a = new FromFileDeconvolutionParameters(feats, 1, 60);
            var b = new FromFileDeconvolutionParameters(feats, 1, 60);
            Assert.That(a, Is.EqualTo(b));
        }

        [Test]
        public void FromFile_Equal_DifferentFeatureCount_AreEqualWhenNoFilePath()
        {
            // Two instances built with the in-memory feature ctor (no FilePath) must
            // be equal under the new config-only equality contract: FilePath is the
            // sole FromFile-specific identity, and feature count is now load-derived
            // state, not configuration. Feature count no longer participates in
            // equality (it would force eager I/O on every hash/equals call).
            var a = new FromFileDeconvolutionParameters(new[] { Feature(), Feature(601.0, 3) }, 1, 60);
            var b = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            Assert.That(a, Is.EqualTo(b));
        }

        [Test]
        public void FromFile_Equal_DifferentFilePaths_AreNotEqual()
        {
            // The new config-only identity is FilePath. Two in-memory instances with
            // no FilePath are equal; assigning different FilePaths must distinguish
            // them. Build in-memory, set FilePath post-construction, then compare.
            var a = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            var b = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            a.FilePath = "featureA.ms1.feature";
            b.FilePath = "featureB.ms1.feature";
            Assert.That(a, Is.Not.EqualTo(b));
        }

        [Test]
        public void FromFile_Equal_SameReference_IsEqual()
        {
            var a = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            Assert.That(a, Is.EqualTo(a));
        }

        [Test]
        public void FromFile_Equal_Null_IsNotEqual()
        {
            var a = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            Assert.That(a.Equals(null), Is.False);
        }

        [Test]
        public void FromFile_Equal_CrossType_IsNotEqual()
        {
            var fromFile = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            var classic = new ClassicDeconvolutionParameters(1, 12, 4.0, 3.0);
            Assert.That(fromFile, Is.Not.EqualTo(classic));
        }

        [Test]
        public void FromFile_Clone_ReturnsNewReference()
        {
            var original = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60);
            var clone = original.Clone();
            Assert.That(clone, Is.Not.SameAs(original));
        }

        [Test]
        public void FromFile_Clone_EqualsOriginal()
        {
            var feats = new[] { Feature(), Feature(601.0, 3, intensity: 2e5) };
            var original = new FromFileDeconvolutionParameters(feats, 1, 60);
            var clone = original.Clone();
            Assert.That(clone, Is.EqualTo(original));
        }

        [Test]
        public void FromFile_Clone_PreservesUseGenericScore()
        {
            var original = new FromFileDeconvolutionParameters(new[] { Feature() }, 1, 60)
            {
                UseGenericScore = true
            };
            var clone = original.Clone();
            Assert.That(clone.UseGenericScore, Is.True);
        }
    }
}
