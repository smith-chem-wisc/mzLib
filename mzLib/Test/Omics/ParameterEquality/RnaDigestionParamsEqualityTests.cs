using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using System;
using Transcriptomics.Digestion;

namespace Test.Omics.ParameterEquality
{
    [TestFixture]
    public class RnaDigestionParamsEqualityTests
    {
        private static RnaDigestionParams Create(
            string rnase = "top-down",
            int maxMissedCleavages = 0,
            int minLength = 3,
            int maxLength = 50,
            int maxModificationIsoforms = 1024,
            int maxMods = 2,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both)
        {
            return new RnaDigestionParams(
                rnase,
                maxMissedCleavages,
                minLength,
                maxLength,
                maxModificationIsoforms,
                maxMods,
                fragmentationTerminus);
        }

        #region RnaDigestionParams Concrete Equality

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var rdp = Create();
            Assert.That(rdp.Equals(rdp), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var rdp = Create();
            Assert.That(rdp.Equals(null!), Is.False);
        }

        [Test]
        public void Equals_IdenticalParams_ReturnsTrue()
        {
            var rdp1 = Create();
            var rdp2 = Create();
            Assert.That(rdp1.Equals(rdp2), Is.True);
        }

        [Test]
        public void Equals_DifferentRnase_ReturnsFalse()
        {
            var rdp1 = Create(rnase: "top-down");
            var rdp2 = Create(rnase: "RNase T1");
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxMissedCleavages_ReturnsFalse()
        {
            var rdp1 = Create(maxMissedCleavages: 0);
            var rdp2 = Create(maxMissedCleavages: 1);
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMinLength_ReturnsFalse()
        {
            var rdp1 = Create(minLength: 3);
            var rdp2 = Create(minLength: 5);
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxLength_ReturnsFalse()
        {
            var rdp1 = Create(maxLength: 50);
            var rdp2 = Create(maxLength: 100);
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxModificationIsoforms_ReturnsFalse()
        {
            var rdp1 = Create(maxModificationIsoforms: 1024);
            var rdp2 = Create(maxModificationIsoforms: 512);
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxMods_ReturnsFalse()
        {
            var rdp1 = Create(maxMods: 2);
            var rdp2 = Create(maxMods: 3);
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        [Test]
        public void Equals_DifferentFragmentationTerminus_ReturnsFalse()
        {
            var rdp1 = Create(fragmentationTerminus: FragmentationTerminus.Both);
            var rdp2 = Create(fragmentationTerminus: FragmentationTerminus.N);
            Assert.That(rdp1.Equals(rdp2), Is.False);
        }

        #endregion

        #region RnaDigestionParams Object Equals

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            var rdp = Create();
            Assert.That(rdp.Equals((object)null), Is.False);
        }

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            var rdp = Create();
            Assert.That(rdp.Equals((object)rdp), Is.True);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            var rdp = Create();
            Assert.That(rdp.Equals("not an rnadigestionparams"), Is.False);
        }

        [Test]
        public void Equals_Object_IdenticalParams_ReturnsTrue()
        {
            var rdp1 = Create();
            var rdp2 = Create();
            Assert.That(rdp1.Equals((object)rdp2), Is.True);
        }

        #endregion

        #region IDigestionParams Interface Equality

        [Test]
        public void Equals_IDigestionParams_SameReference_ReturnsTrue()
        {
            var rdp = Create();
            IDigestionParams idp = rdp;
            Assert.That(rdp.Equals(idp), Is.True);
        }

        [Test]
        public void Equals_IDigestionParams_Identical_ReturnsTrue()
        {
            var rdp1 = Create();
            var rdp2 = Create();
            IDigestionParams idp1 = rdp1;
            IDigestionParams idp2 = rdp2;
            Assert.That(idp1.Equals(idp2), Is.True);
        }

        [Test]
        public void Equals_IDigestionParams_DifferentMaxMissedCleavages_ReturnsFalse()
        {
            var rdp1 = Create(maxMissedCleavages: 0);
            var rdp2 = Create(maxMissedCleavages: 1);
            IDigestionParams idp1 = rdp1;
            IDigestionParams idp2 = rdp2;
            Assert.That(idp1.Equals(idp2), Is.False);
        }

        [Test]
        public void Equals_IDigestionParams_CrossType_ReturnsFalse()
        {
            var rdp = Create();
            IDigestionParams idp = new Proteomics.ProteolyticDigestion.DigestionParams();
            Assert.That(rdp.Equals(idp), Is.False);
        }

        #endregion

        #region RnaDigestionParams GetHashCode

        [Test]
        public void GetHashCode_IdenticalParams_AreEqual()
        {
            var rdp1 = Create();
            var rdp2 = Create();
            Assert.That(rdp1.GetHashCode(), Is.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentRnase_AreNotEqual()
        {
            var rdp1 = Create(rnase: "top-down");
            var rdp2 = Create(rnase: "RNase T1");
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentMaxMissedCleavages_AreNotEqual()
        {
            var rdp1 = Create(maxMissedCleavages: 0);
            var rdp2 = Create(maxMissedCleavages: 1);
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentMinLength_AreNotEqual()
        {
            var rdp1 = Create(minLength: 3);
            var rdp2 = Create(minLength: 5);
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentMaxLength_AreNotEqual()
        {
            var rdp1 = Create(maxLength: 50);
            var rdp2 = Create(maxLength: 100);
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentMaxModificationIsoforms_AreNotEqual()
        {
            var rdp1 = Create(maxModificationIsoforms: 1024);
            var rdp2 = Create(maxModificationIsoforms: 512);
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentMaxMods_AreNotEqual()
        {
            var rdp1 = Create(maxMods: 2);
            var rdp2 = Create(maxMods: 3);
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFragmentationTerminus_AreNotEqual()
        {
            var rdp1 = Create(fragmentationTerminus: FragmentationTerminus.Both);
            var rdp2 = Create(fragmentationTerminus: FragmentationTerminus.N);
            Assert.That(rdp1.GetHashCode(), Is.Not.EqualTo(rdp2.GetHashCode()));
        }

        #endregion
    }
}