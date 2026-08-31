using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using Transcriptomics.Digestion;

namespace Test.Omics.ParameterEquality
{
    [TestFixture]
    public class DigestionParamsEqualityTests
    {
        private static DigestionParams Create(
            string protease = "trypsin",
            int maxMissedCleavages = 2,
            int minPeptideLength = 7,
            int maxPeptideLength = 50,
            int maxModificationIsoforms = 1024,
            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
            int maxModsForPeptides = 2,
            CleavageSpecificity searchModeType = CleavageSpecificity.Full,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            bool generateUnlabeledProteinsForSilac = true,
            bool keepNGlycopeptide = false,
            bool keepOGlycopeptide = false)
        {
            return new DigestionParams(
                protease,
                maxMissedCleavages,
                minPeptideLength,
                maxPeptideLength,
                maxModificationIsoforms,
                initiatorMethionineBehavior,
                maxModsForPeptides,
                searchModeType,
                fragmentationTerminus,
                generateUnlabeledProteinsForSilac,
                keepNGlycopeptide,
                keepOGlycopeptide);
        }

        #region DigestionParams Concrete Equality

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var dp = Create();
            Assert.That(dp.Equals(dp), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var dp = Create();
            Assert.That(dp.Equals((DigestionParams)null!), Is.False);
        }

        [Test]
        public void Equals_IdenticalParams_ReturnsTrue()
        {
            var dp1 = Create();
            var dp2 = Create();
            Assert.That(dp1.Equals(dp2), Is.True);
        }

        [Test]
        public void Equals_DifferentProtease_ReturnsFalse()
        {
            var dp1 = Create(protease: "trypsin");
            var dp2 = Create(protease: "Arg-C");
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxMissedCleavages_ReturnsFalse()
        {
            var dp1 = Create(maxMissedCleavages: 2);
            var dp2 = Create(maxMissedCleavages: 3);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMinPeptideLength_ReturnsFalse()
        {
            var dp1 = Create(minPeptideLength: 7);
            var dp2 = Create(minPeptideLength: 8);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxPeptideLength_ReturnsFalse()
        {
            var dp1 = Create(maxPeptideLength: 50);
            var dp2 = Create(maxPeptideLength: 100);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentInitiatorMethionineBehavior_ReturnsFalse()
        {
            var dp1 = Create(initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable);
            var dp2 = Create(initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxModificationIsoforms_ReturnsFalse()
        {
            var dp1 = Create(maxModificationIsoforms: 1024);
            var dp2 = Create(maxModificationIsoforms: 512);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMaxMods_ReturnsFalse()
        {
            var dp1 = Create(maxModsForPeptides: 2);
            var dp2 = Create(maxModsForPeptides: 3);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentSearchModeType_ReturnsFalse()
        {
            var dp1 = Create(searchModeType: CleavageSpecificity.Full);
            var dp2 = Create(searchModeType: CleavageSpecificity.Semi);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentFragmentationTerminus_ReturnsFalse()
        {
            var dp1 = Create(fragmentationTerminus: FragmentationTerminus.Both);
            var dp2 = Create(fragmentationTerminus: FragmentationTerminus.N);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentGenerateUnlabeledProteinsForSilac_ReturnsFalse()
        {
            var dp1 = Create(generateUnlabeledProteinsForSilac: true);
            var dp2 = Create(generateUnlabeledProteinsForSilac: false);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentKeepNGlycopeptide_ReturnsFalse()
        {
            var dp1 = Create(keepNGlycopeptide: false);
            var dp2 = Create(keepNGlycopeptide: true);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        [Test]
        public void Equals_DifferentKeepOGlycopeptide_ReturnsFalse()
        {
            var dp1 = Create(keepOGlycopeptide: false);
            var dp2 = Create(keepOGlycopeptide: true);
            Assert.That(dp1.Equals(dp2), Is.False);
        }

        #endregion

        #region DigestionParams Object Equals

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            var dp = Create();
            Assert.That(dp.Equals((object?)null), Is.False);
        }

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            var dp = Create();
            Assert.That(dp.Equals((object)dp), Is.True);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            var dp = Create();
            Assert.That(dp.Equals("not a digestionparams"), Is.False);
        }

        [Test]
        public void Equals_Object_IdenticalParams_ReturnsTrue()
        {
            var dp1 = Create();
            var dp2 = Create();
            Assert.That(dp1.Equals((object)dp2), Is.True);
        }

        #endregion

        #region IDigestionParams Interface Equality

        [Test]
        public void Equals_IDigestionParams_SameReference_ReturnsTrue()
        {
            var dp = Create();
            IDigestionParams idp = dp;
            Assert.That(dp.Equals(idp), Is.True);
        }

        [Test]
        public void Equals_IDigestionParams_Identical_ReturnsTrue()
        {
            var dp1 = Create();
            var dp2 = Create();
            IDigestionParams idp1 = dp1;
            IDigestionParams idp2 = dp2;
            Assert.That(idp1.Equals(idp2), Is.True);
        }

        [Test]
        public void Equals_IDigestionParams_DifferentMaxMissedCleavages_ReturnsFalse()
        {
            var dp1 = Create(maxMissedCleavages: 2);
            var dp2 = Create(maxMissedCleavages: 3);
            IDigestionParams idp1 = dp1;
            IDigestionParams idp2 = dp2;
            Assert.That(idp1.Equals(idp2), Is.False);
        }

        [Test]
        public void Equals_IDigestionParams_CrossType_ReturnsFalse()
        {
            var dp = Create();
            IDigestionParams idp = new RnaDigestionParams();
            Assert.That(dp.Equals(idp), Is.False);
        }

        #endregion

        #region DigestionParams GetHashCode

        [Test]
        public void GetHashCode_IdenticalParams_AreEqual()
        {
            var dp1 = Create();
            var dp2 = Create();
            Assert.That(dp1.GetHashCode(), Is.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentProtease_AreNotEqual()
        {
            var dp1 = Create(protease: "trypsin");
            var dp2 = Create(protease: "Arg-C");
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentMaxMissedCleavages_AreNotEqual()
        {
            var dp1 = Create(maxMissedCleavages: 2);
            var dp2 = Create(maxMissedCleavages: 3);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentInitiatorMethionineBehavior_AreNotEqual()
        {
            var dp1 = Create(initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable);
            var dp2 = Create(initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentSearchModeType_AreNotEqual()
        {
            var dp1 = Create(searchModeType: CleavageSpecificity.Full);
            var dp2 = Create(searchModeType: CleavageSpecificity.Semi);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentFragmentationTerminus_AreNotEqual()
        {
            var dp1 = Create(fragmentationTerminus: FragmentationTerminus.Both);
            var dp2 = Create(fragmentationTerminus: FragmentationTerminus.N);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentGenerateUnlabeledProteinsForSilac_AreNotEqual()
        {
            var dp1 = Create(generateUnlabeledProteinsForSilac: true);
            var dp2 = Create(generateUnlabeledProteinsForSilac: false);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentKeepNGlycopeptide_AreNotEqual()
        {
            var dp1 = Create(keepNGlycopeptide: false);
            var dp2 = Create(keepNGlycopeptide: true);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentKeepOGlycopeptide_AreNotEqual()
        {
            var dp1 = Create(keepOGlycopeptide: false);
            var dp2 = Create(keepOGlycopeptide: true);
            Assert.That(dp1.GetHashCode(), Is.Not.EqualTo(dp2.GetHashCode()));
        }

        #endregion
    }
}