using NUnit.Framework;
using Omics.Fragmentation;
using System.Collections.Generic;
using Transcriptomics;
using FragmentationParams = Omics.Fragmentation.FragmentationParams;

namespace Test.Omics.ParameterEquality
{
    [TestFixture]
    public class RnaFragmentationParamsEqualityTests
    {
        private static RnaFragmentationParams Create(
            bool generateMIon = true,
            List<MIonLoss> mIonLosses = null,
            bool modificationsCanSuppressBaseLossIons = false)
        {
            return new RnaFragmentationParams
            {
                GenerateMIon = generateMIon,
                MIonLosses = mIonLosses ?? new List<MIonLoss>(),
                ModificationsCanSuppressBaseLossIons = modificationsCanSuppressBaseLossIons
            };
        }

        #region RnaFragmentationParams Concrete Equality

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var rfp = Create();
            Assert.That(rfp.Equals(rfp), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var rfp = Create();
            Assert.That(rfp.Equals(null!), Is.False);
        }

        [Test]
        public void Equals_IdenticalParams_ReturnsTrue()
        {
            var rfp1 = Create();
            var rfp2 = Create();
            Assert.That(rfp1.Equals(rfp2), Is.True);
        }

        [Test]
        public void Equals_DifferentGenerateMIon_ReturnsFalse()
        {
            var rfp1 = Create(generateMIon: true);
            var rfp2 = Create(generateMIon: false);
            Assert.That(rfp1.Equals(rfp2), Is.False);
        }

        [Test]
        public void Equals_DifferentModificationsCanSuppressBaseLossIons_ReturnsFalse()
        {
            var rfp1 = Create(modificationsCanSuppressBaseLossIons: true);
            var rfp2 = Create(modificationsCanSuppressBaseLossIons: false);
            Assert.That(rfp1.Equals(rfp2), Is.False);
        }

        [Test]
        public void Equals_DifferentMIonLosses_ReturnsFalse()
        {
            var rfp1 = Create(mIonLosses: new List<MIonLoss> { new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O")) });
            var rfp2 = Create(mIonLosses: new List<MIonLoss> { new MIonLoss("B", "-B", Chemistry.ChemicalFormula.ParseFormula("CO2")) });
            Assert.That(rfp1.Equals(rfp2), Is.False);
        }

        [Test]
        public void Equals_IdenticalMIonLosses_ReturnsTrue()
        {
            var loss1 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var loss2 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var rfp1 = Create(mIonLosses: new List<MIonLoss> { loss1 });
            var rfp2 = Create(mIonLosses: new List<MIonLoss> { loss2 });
            Assert.That(rfp1.Equals(rfp2), Is.True);
        }

        [Test]
        public void Equals_MIonLossesDifferentOrder_ReturnsTrue()
        {
            var lossA = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var lossB = new MIonLoss("B", "-B", Chemistry.ChemicalFormula.ParseFormula("CO2"));
            var rfp1 = Create(mIonLosses: new List<MIonLoss> { lossA, lossB });
            var rfp2 = Create(mIonLosses: new List<MIonLoss> { lossB, lossA });
            Assert.That(rfp1.Equals(rfp2), Is.True);
        }

        #endregion

        #region RnaFragmentationParams Object Equals

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            var rfp = Create();
            Assert.That(rfp.Equals((object)null), Is.False);
        }

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            var rfp = Create();
            Assert.That(rfp.Equals((object)rfp), Is.True);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            var rfp = Create();
            Assert.That(rfp.Equals("not an rnafragmentationparams"), Is.False);
        }

        [Test]
        public void Equals_Object_IdenticalParams_ReturnsTrue()
        {
            var rfp1 = Create();
            var rfp2 = Create();
            Assert.That(rfp1.Equals((object)rfp2), Is.True);
        }

        #endregion

        #region IFragmentationParams Interface Equality

        [Test]
        public void Equals_IFragmentationParams_SameReference_ReturnsTrue()
        {
            var rfp = Create();
            IFragmentationParams ifp = rfp;
            Assert.That(rfp.Equals(ifp), Is.True);
        }

        [Test]
        public void Equals_IFragmentationParams_Identical_ReturnsTrue()
        {
            var rfp1 = Create();
            var rfp2 = Create();
            IFragmentationParams ifp1 = rfp1;
            IFragmentationParams ifp2 = rfp2;
            Assert.That(ifp1.Equals(ifp2), Is.True);
        }

        [Test]
        public void Equals_IFragmentationParams_DifferentGenerateMIon_ReturnsFalse()
        {
            var rfp1 = Create(generateMIon: true);
            var rfp2 = Create(generateMIon: false);
            IFragmentationParams ifp1 = rfp1;
            IFragmentationParams ifp2 = rfp2;
            Assert.That(ifp1.Equals(ifp2), Is.False);
        }

        [Test]
        public void Equals_IFragmentationParams_CrossType_ReturnsFalse()
        {
            var rfp = Create();
            IFragmentationParams ifp = new FragmentationParams();
            Assert.That(rfp.Equals(ifp), Is.False);
        }

        [Test]
        public void Equals_IFragmentationParams_ViaInterface_IncludesModificationsCanSuppressBaseLossIons()
        {
            var rfp1 = Create(modificationsCanSuppressBaseLossIons: true);
            var rfp2 = Create(modificationsCanSuppressBaseLossIons: false);
            IFragmentationParams ifp1 = rfp1;
            IFragmentationParams ifp2 = rfp2;
            Assert.That(ifp1.Equals(ifp2), Is.False);
        }

        #endregion

        #region RnaFragmentationParams GetHashCode

        [Test]
        public void GetHashCode_IdenticalParams_AreEqual()
        {
            var rfp1 = Create();
            var rfp2 = Create();
            Assert.That(rfp1.GetHashCode(), Is.EqualTo(rfp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentParams_AreNotEqual()
        {
            var rfp1 = Create(generateMIon: true);
            var rfp2 = Create(generateMIon: false);
            Assert.That(rfp1.GetHashCode(), Is.Not.EqualTo(rfp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentModificationsCanSuppressBaseLossIons_AreNotEqual()
        {
            var rfp1 = Create(modificationsCanSuppressBaseLossIons: true);
            var rfp2 = Create(modificationsCanSuppressBaseLossIons: false);
            Assert.That(rfp1.GetHashCode(), Is.Not.EqualTo(rfp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_IdenticalMIonLosses_AreEqual()
        {
            var loss1 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var loss2 = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var rfp1 = Create(mIonLosses: new List<MIonLoss> { loss1 });
            var rfp2 = Create(mIonLosses: new List<MIonLoss> { loss2 });
            Assert.That(rfp1.GetHashCode(), Is.EqualTo(rfp2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_MIonLossesDifferentOrder_AreEqual()
        {
            var lossA = new MIonLoss("A", "-A", Chemistry.ChemicalFormula.ParseFormula("H2O"));
            var lossB = new MIonLoss("B", "-B", Chemistry.ChemicalFormula.ParseFormula("CO2"));
            var rfp1 = Create(mIonLosses: new List<MIonLoss> { lossA, lossB });
            var rfp2 = Create(mIonLosses: new List<MIonLoss> { lossB, lossA });
            Assert.That(rfp1.GetHashCode(), Is.EqualTo(rfp2.GetHashCode()));
        }

        #endregion
    }
}