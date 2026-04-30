using NUnit.Framework;
using Omics.Fragmentation;
using System.Collections.Generic;

namespace Test.Omics.ParameterEquality
{
    [TestFixture]
    public class FragmentationParamsEqualityTests
    {
        private static FragmentationParams Create(bool generateMIon = false, List<MIonLoss>? mIonLosses = null)
        {
            return new FragmentationParams
            {
                GenerateMIon = generateMIon,
                MIonLosses = mIonLosses ?? new List<MIonLoss>()
            };
        }

        #region FragmentationParams Concrete Equality

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var fp = Create();
            Assert.That(fp.Equals(fp), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var fp = Create();
            Assert.That(fp.Equals((FragmentationParams)null!), Is.False);
        }

        [Test]
        public void Equals_IdenticalParams_ReturnsTrue()
        {
            var fp1 = Create(generateMIon: false);
            var fp2 = Create(generateMIon: false);
            Assert.That(fp1.Equals(fp2), Is.True);
        }

        [Test]
        public void Equals_DifferentGenerateMIon_ReturnsFalse()
        {
            var fp1 = Create(generateMIon: true);
            var fp2 = Create(generateMIon: false);
            Assert.That(fp1.Equals(fp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMIonLosses_ReturnsFalse()
        {
            var fp1 = Create(mIonLosses: new List<MIonLoss> { new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O")) });
            var fp2 = Create(mIonLosses: new List<MIonLoss> { new MIonLoss("B", "-B", Chemistry.ChemicalFormula.ParseFormula("CO2")) });
            Assert.That(fp1.Equals(fp2), Is.False);
        }

        [Test]
        public void Equals_IdenticalMIonLosses_ReturnsTrue()
        {
            var loss1 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var loss2 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var fp1 = Create(mIonLosses: new List<MIonLoss> { loss1 });
            var fp2 = Create(mIonLosses: new List<MIonLoss> { loss2 });
            Assert.That(fp1.Equals(fp2), Is.True);
        }

        [Test]
        public void Equals_MIonLossesDifferentOrder_ReturnsTrue()
        {
            var lossA = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var lossB = new MIonLoss("B", "-B", Chemistry.ChemicalFormula.ParseFormula("CO2"));
            var fp1 = Create(mIonLosses: new List<MIonLoss> { lossA, lossB });
            var fp2 = Create(mIonLosses: new List<MIonLoss> { lossB, lossA });
            Assert.That(fp1.Equals(fp2), Is.True);
        }

        #endregion

        #region FragmentationParams Object Equals

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            var fp = Create();
            Assert.That(fp.Equals((object?)null), Is.False);
        }

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            var fp = Create();
            Assert.That(fp.Equals((object)fp), Is.True);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            var fp = Create();
            Assert.That(fp.Equals("not a fragmentationparams"), Is.False);
        }

        [Test]
        public void Equals_Object_IdenticalParams_ReturnsTrue()
        {
            var fp1 = Create(generateMIon: false);
            var fp2 = Create(generateMIon: false);
            Assert.That(fp1.Equals((object)fp2), Is.True);
        }

        #endregion

        #region IFragmentationParams Interface Equality

        [Test]
        public void Equals_IFragmentationParams_SameReference_ReturnsTrue()
        {
            var fp = Create();
            IFragmentationParams ifp = fp;
            Assert.That(fp.Equals(ifp), Is.True);
        }

        [Test]
        public void Equals_IFragmentationParams_Identical_ReturnsTrue()
        {
            var fp1 = Create(generateMIon: false);
            var fp2 = Create(generateMIon: false);
            IFragmentationParams ifp1 = fp1;
            IFragmentationParams ifp2 = fp2;
            Assert.That(ifp1.Equals(ifp2), Is.True);
        }

        [Test]
        public void Equals_IFragmentationParams_DifferentGenerateMIon_ReturnsFalse()
        {
            var fp1 = Create(generateMIon: true);
            var fp2 = Create(generateMIon: false);
            IFragmentationParams ifp1 = fp1;
            IFragmentationParams ifp2 = fp2;
            Assert.That(ifp1.Equals(ifp2), Is.False);
        }



        #endregion

        #region FragmentationParams GetHashCode

        [Test]
        public void GetHashCode_IdenticalParams_AreEqual()
        {
            var fp1 = Create(generateMIon: false);
            var fp2 = Create(generateMIon: false);
            Assert.That(fp1.GetHashCode(), Is.EqualTo(fp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentParams_AreNotEqual()
        {
            var fp1 = Create(generateMIon: true);
            var fp2 = Create(generateMIon: false);
            Assert.That(fp1.GetHashCode(), Is.Not.EqualTo(fp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_IdenticalMIonLosses_AreEqual()
        {
            var loss1 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var loss2 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var fp1 = Create(mIonLosses: new List<MIonLoss> { loss1 });
            var fp2 = Create(mIonLosses: new List<MIonLoss> { loss2 });
            Assert.That(fp1.GetHashCode(), Is.EqualTo(fp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_MIonLossesDifferentOrder_AreEqual()
        {
            var lossA = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var lossB = new MIonLoss("B", "-B", Chemistry.ChemicalFormula.ParseFormula("CO2"));
            var fp1 = Create(mIonLosses: new List<MIonLoss> { lossA, lossB });
            var fp2 = Create(mIonLosses: new List<MIonLoss> { lossB, lossA });
            Assert.That(fp1.GetHashCode(), Is.EqualTo(fp2.GetHashCode()));
        }

        #endregion
    }
}