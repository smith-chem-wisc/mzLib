using Chemistry;
using NUnit.Framework;
using Omics.Fragmentation;
using System.Collections.Generic;

namespace Test.Omics.ParameterEquality
{
    [TestFixture]
    public class MIonLossEqualityTests
    {
        private static MIonLoss Create(string name, string annotation, string formula)
            => new MIonLoss(name, annotation, ChemicalFormula.ParseFormula(formula));

        #region MIonLoss Equality

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var loss = Create("Test", "-T", "H2O");
            Assert.That(loss.Equals(loss), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var loss = Create("Test", "-T", "H2O");
            Assert.That(loss.Equals((MIonLoss)null!), Is.False);
        }

        [Test]
        public void Equals_DifferentName_ReturnsFalse()
        {
            var loss1 = Create("Name1", "-T", "H2O");
            var loss2 = Create("Name2", "-T", "H2O");
            Assert.That(loss1.Equals(loss2), Is.False);
        }

        [Test]
        public void Equals_DifferentAnnotation_ReturnsFalse()
        {
            var loss1 = Create("Name", "-A", "H2O");
            var loss2 = Create("Name", "-B", "H2O");
            Assert.That(loss1.Equals(loss2), Is.False);
        }

        [Test]
        public void Equals_DifferentChemicalFormula_ReturnsFalse()
        {
            var loss1 = Create("Name", "-T", "H2O");
            var loss2 = Create("Name", "-T", "CO2");
            Assert.That(loss1.Equals(loss2), Is.False);
        }

        [Test]
        public void Equals_IdenticalValues_ReturnsTrue()
        {
            var loss1 = Create("Name", "-T", "H2O");
            var loss2 = Create("Name", "-T", "H2O");
            Assert.That(loss1.Equals(loss2), Is.True);
        }

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            var loss = Create("Test", "-T", "H2O");
            Assert.That(loss.Equals((object?)null), Is.False);
        }

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            var loss = Create("Test", "-T", "H2O");
            Assert.That(loss.Equals((object)loss), Is.True);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            var loss = Create("Test", "-T", "H2O");
            Assert.That(loss.Equals("not an mionloss"), Is.False);
        }

        [Test]
        public void Equals_Object_IdenticalValues_ReturnsTrue()
        {
            var loss1 = Create("Name", "-T", "H2O");
            var loss2 = Create("Name", "-T", "H2O");
            Assert.That(loss1.Equals((object)loss2), Is.True);
        }

        [Test]
        public void GetHashCode_IdenticalObjects_AreEqual()
        {
            var loss1 = Create("Name", "-T", "H2O");
            var loss2 = Create("Name", "-T", "H2O");
            Assert.That(loss1.GetHashCode(), Is.EqualTo(loss2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentObjects_AreNotEqual()
        {
            var loss1 = Create("Name", "-T", "H2O");
            var loss2 = Create("Name", "-T", "CO2");
            Assert.That(loss1.GetHashCode(), Is.Not.EqualTo(loss2.GetHashCode()));
        }

        #endregion

        #region MIonListComparer

        [Test]
        public void MIonListComparer_Equals_SameOrder_ReturnsTrue()
        {
            var list1 = new List<MIonLoss> { Create("A", "-A", "H2O"), Create("B", "-B", "CO2") };
            var list2 = new List<MIonLoss> { Create("A", "-A", "H2O"), Create("B", "-B", "CO2") };
            Assert.That(MIonListComparer.Instance.Equals(list1, list2), Is.True);
        }

        [Test]
        public void MIonListComparer_Equals_DifferentOrder_ReturnsTrue()
        {
            var list1 = new List<MIonLoss> { Create("A", "-A", "H2O"), Create("B", "-B", "CO2") };
            var list2 = new List<MIonLoss> { Create("B", "-B", "CO2"), Create("A", "-A", "H2O") };
            Assert.That(MIonListComparer.Instance.Equals(list1, list2), Is.True);
        }

        [Test]
        public void MIonListComparer_Equals_DifferentCount_ReturnsFalse()
        {
            var list1 = new List<MIonLoss> { Create("A", "-A", "H2O"), Create("B", "-B", "CO2") };
            var list2 = new List<MIonLoss> { Create("A", "-A", "H2O") };
            Assert.That(MIonListComparer.Instance.Equals(list1, list2), Is.False);
        }

        [Test]
        public void MIonListComparer_Equals_DifferentItems_ReturnsFalse()
        {
            var list1 = new List<MIonLoss> { Create("A", "-A", "H2O") };
            var list2 = new List<MIonLoss> { Create("B", "-B", "CO2") };
            Assert.That(MIonListComparer.Instance.Equals(list1, list2), Is.False);
        }

        [Test]
        public void MIonListComparer_Equals_NullX_ReturnsFalse()
        {
            var list = new List<MIonLoss> { Create("A", "-A", "H2O") };
            Assert.That(MIonListComparer.Instance.Equals(null!, list), Is.False);
        }

        [Test]
        public void MIonListComparer_Equals_NullY_ReturnsFalse()
        {
            var list = new List<MIonLoss> { Create("A", "-A", "H2O") };
            Assert.That(MIonListComparer.Instance.Equals(list, null!), Is.False);
        }

        [Test]
        public void MIonListComparer_Equals_BothNull_ReturnsFalse()
        {
            Assert.That(MIonListComparer.Instance.Equals(null!, null!), Is.False);
        }

        [Test]
        public void MIonListComparer_GetHashCode_SameListDifferentOrder_AreEqual()
        {
            var list1 = new List<MIonLoss> { Create("A", "-A", "H2O"), Create("B", "-B", "CO2") };
            var list2 = new List<MIonLoss> { Create("B", "-B", "CO2"), Create("A", "-A", "H2O") };
            Assert.That(MIonListComparer.Instance.GetHashCode(list1), Is.EqualTo(MIonListComparer.Instance.GetHashCode(list2)));
        }

        [Test]
        public void MIonListComparer_GetHashCode_DifferentLists_AreNotEqual()
        {
            var list1 = new List<MIonLoss> { Create("A", "-A", "H2O") };
            var list2 = new List<MIonLoss> { Create("B", "-B", "CO2") };
            Assert.That(MIonListComparer.Instance.GetHashCode(list1), Is.Not.EqualTo(MIonListComparer.Instance.GetHashCode(list2)));
        }

        [Test]
        public void MIonListComparer_Instance_IsSingleton()
        {
            Assert.That(MIonListComparer.Instance, Is.SameAs(MIonListComparer.Instance));
        }

        #endregion
    }
}