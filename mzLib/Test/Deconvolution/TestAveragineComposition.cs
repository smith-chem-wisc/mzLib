using MassSpectrometry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Collections.Generic;

namespace Test.Deconvolution
{
    /// <summary>
    /// Tests for the public average-composition accessors on <see cref="Averagine"/> and
    /// <see cref="OxyriboAveragine"/>, which expose the per-residue elemental composition so callers
    /// (e.g. FlashLFQ) can derive an approximate chemical formula for an arbitrary mass.
    /// </summary>
    [TestFixture]
    public static class TestAveragineComposition
    {
        [Test]
        public static void AveragineCompositionMatchesSenkoValues()
        {
            var composition = new Averagine().GetAverageChemicalFormula();

            Assert.AreEqual(5, composition.Count);
            Assert.AreEqual(4.9384, composition['C'], 1e-9);
            Assert.AreEqual(7.7583, composition['H'], 1e-9);
            Assert.AreEqual(1.4773, composition['O'], 1e-9);
            Assert.AreEqual(1.3577, composition['N'], 1e-9);
            Assert.AreEqual(0.0417, composition['S'], 1e-9);
        }

        [Test]
        public static void OxyriboAveragineCompositionIsMeanRibonucleotide()
        {
            var composition = new OxyriboAveragine().GetAverageChemicalFormula();

            // mean of AMP/CMP/GMP/UMP: C9.5 H13.75 N3.75 O8 P1
            Assert.AreEqual(5, composition.Count);
            Assert.AreEqual(9.5, composition['C'], 1e-9);
            Assert.AreEqual(13.75, composition['H'], 1e-9);
            Assert.AreEqual(8.0, composition['O'], 1e-9);
            Assert.AreEqual(3.75, composition['N'], 1e-9);
            Assert.AreEqual(1.0, composition['P'], 1e-9);
        }

        [Test]
        public static void GetAverageChemicalFormulaReturnsIndependentCopy()
        {
            // mutating a returned dictionary must not corrupt the shared static composition
            var first = new Averagine().GetAverageChemicalFormula();
            first['C'] = -1;

            var second = new Averagine().GetAverageChemicalFormula();
            Assert.AreEqual(4.9384, second['C'], 1e-9);
        }
    }
}
