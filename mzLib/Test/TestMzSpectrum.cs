using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownProteomics.IO.MzIdentMl;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzSpectrum
    {

        [Test]
        public void Bubba()
        {
            double[] mzs = new double[3] { 1,1.0000001, 1.0000002 };
            double[] intensities = new double[3] { 1,2,3 };
            bool shouldCopy = false;
            MzSpectrum s = new MzSpectrum(mzs,intensities,shouldCopy);

            s.MergeMzAndIntensityValuesForAnyMzsWithinTolerance(new PpmTolerance(0.01));
        }

        [Test]
        public void Bubba2()
        {
            var j = ExtensionMethods.MergeTwoSortedArraysIntoSortedArray(new double[] { 2, 4, 6, }, new double[] { 1, 3, 5, 7, 8, 9, 10, 11, 12 });
            CollectionAssert.AreEquivalent(new double[]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }, j);
        }
    }
}
