using NUnit.Framework;
using MassSpectrometry;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public class IsotopicEnvelopeEqualityTests
    {
        private IsotopicEnvelope CreateEnvelope(
            List<(double mz, double intensity)> peaks = null,
            double monoMass = 100.0,
            int charge = 2,
            double totalIntensity = 200.0,
            double stDev = 1.0)
        {
            peaks ??= new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) };
            return new IsotopicEnvelope(peaks, monoMass, charge, totalIntensity, stDev);
        }

        [Test]
        public void Equals_SameReference_ReturnsTrue()
        {
            var env = CreateEnvelope();
            Assert.That(env.Equals(env), Is.True);
        }

        [Test]
        public void Equals_Null_ReturnsFalse()
        {
            var env = CreateEnvelope();
            Assert.That(env.Equals((IsotopicEnvelope)null), Is.False);
        }

        [Test]
        public void Equals_DifferentCharge_ReturnsFalse()
        {
            var env1 = CreateEnvelope(charge: 2);
            var env2 = CreateEnvelope(charge: 3);
            Assert.That(env1.Equals(env2), Is.False);
        }

        [Test]
        public void Equals_DifferentPeaksCount_ReturnsFalse()
        {
            var env1 = CreateEnvelope(peaks: new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) });
            var env2 = CreateEnvelope(peaks: new List<(double mz, double intensity)> { (100.0, 50.0) });
            Assert.That(env1.Equals(env2), Is.False);
        }

        [Test]
        public void Equals_DifferentTotalIntensity_ReturnsFalse()
        {
            var env1 = CreateEnvelope(totalIntensity: 200.0);
            var env2 = CreateEnvelope(totalIntensity: 201.0);
            Assert.That(env1.Equals(env2), Is.False);
        }

        [Test]
        public void Equals_DifferentMonoisotopicMass_ReturnsFalse()
        {
            var env1 = CreateEnvelope(monoMass: 100.0);
            var env2 = CreateEnvelope(monoMass: 101.0);
            Assert.That(env1.Equals(env2), Is.False);
        }

        [Test]
        public void Equals_DifferentMostAbundantIsotopicMass_ReturnsFalse()
        {
            var env1 = CreateEnvelope(peaks: new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) });
            var env2 = CreateEnvelope(peaks: new List<(double mz, double intensity)> { (102.0, 50.0), (103.0, 150.0) });
            Assert.That(env1.Equals(env2), Is.False);
        }

        [Test]
        public void Equals_DifferentPeakValues_ReturnsFalse()
        {
            var env1 = CreateEnvelope(peaks: new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) });
            var env2 = CreateEnvelope(peaks: new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 151.0) });
            Assert.That(env1.Equals(env2), Is.False);
        }

        [Test]
        public void Equals_IdenticalValues_ReturnsTrue()
        {
            var peaks = new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) };
            var env1 = CreateEnvelope(peaks: peaks, monoMass: 100.0, charge: 2, totalIntensity: 200.0);
            var env2 = CreateEnvelope(peaks: new List<(double mz, double intensity)>(peaks), monoMass: 100.0, charge: 2, totalIntensity: 200.0);
            Assert.That(env1.Equals(env2), Is.True);
        }

        [Test]
        public void Equals_Object_Null_ReturnsFalse()
        {
            var env = CreateEnvelope();
            Assert.That(env.Equals((object)null), Is.False);
        }

        [Test]
        public void Equals_Object_SameReference_ReturnsTrue()
        {
            var env = CreateEnvelope();
            Assert.That(env.Equals((object)env), Is.True);
        }

        [Test]
        public void Equals_Object_DifferentType_ReturnsFalse()
        {
            var env = CreateEnvelope();
            Assert.That(env.Equals("not an envelope"), Is.False);
        }

        [Test]
        public void Equals_Object_IdenticalValues_ReturnsTrue()
        {
            var peaks = new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) };
            var env1 = CreateEnvelope(peaks: peaks, monoMass: 100.0, charge: 2, totalIntensity: 200.0);
            var env2 = CreateEnvelope(peaks: new List<(double mz, double intensity)>(peaks), monoMass: 100.0, charge: 2, totalIntensity: 200.0);
            Assert.That(env1.Equals((object)env2), Is.True);
        }

        [Test]
        public void GetHashCode_IdenticalObjects_AreEqual()
        {
            var peaks = new List<(double mz, double intensity)> { (100.0, 50.0), (101.0, 150.0) };
            var env1 = CreateEnvelope(peaks: peaks, monoMass: 100.0, charge: 2, totalIntensity: 200.0);
            var env2 = CreateEnvelope(peaks: new List<(double mz, double intensity)>(peaks), monoMass: 100.0, charge: 2, totalIntensity: 200.0);
            Assert.That(env1.GetHashCode(), Is.EqualTo(env2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentObjects_AreNotEqual()
        {
            var env1 = CreateEnvelope(monoMass: 100.0);
            var env2 = CreateEnvelope(monoMass: 101.0);
            Assert.That(env1.GetHashCode(), Is.Not.EqualTo(env2.GetHashCode()));
        }
    }
}
