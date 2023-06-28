using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using NUnit.Framework;
using Transcriptomics;

namespace Test.Transcriptomics
{
    /// <summary>
    /// Test Data generated with  http://rna.rega.kuleuven.be/masspec/mongo.htm
    /// </summary>
    [TestFixture]
    internal class TestNucleicAcid
    {

        [Test]
        [TestCase("GUACUG", 1954.247)]
        [TestCase("A", 347.062)]
        [TestCase("C", 323.051)]
        [TestCase("U", 324.035)]
        [TestCase("G", 363.057)]
        [TestCase("GU", 669.082)]
        [TestCase("AAA", 1005.166)]
        [TestCase("CCC", 933.133)]
        [TestCase("UUU", 936.085)]
        [TestCase("GGG", 1053.151)]
        public void Test(string sequence, double monoMass)
        {
            RNA rna = new RNA(sequence);

            Assert.IsNotNull(rna);
            Assert.That(rna.Length, Is.EqualTo(sequence.Length));
            Assert.That(rna.MonoisotopicMass, Is.EqualTo(monoMass).Within(0.01));

            var rna2 = new RNA(sequence, NucleicAcid.DefaultFivePrimeTerminus, NucleicAcid.DefaultThreePrimeTerminus);

            Assert.IsNotNull(rna2);
            Assert.That(rna2.Length, Is.EqualTo(sequence.Length));
            Assert.That(rna2.MonoisotopicMass, Is.EqualTo(monoMass).Within(0.01));

            Assert.That(rna.Equals(rna2));
            Assert.That(rna.Equals((object)rna2));
        }
    }
}
