using NUnit.Framework;
using Proteomics.PSM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics.Fragmentation;
using Omics.SpectrumMatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestLibrarySpectrum
    {
        [Test]
        public static void TestDecoyLibrarySpectraGenerationFunction()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            Product c = new Product(ProductType.b, FragmentationTerminus.N, 3, 3, 1, 0);
            Product d = new Product(ProductType.b, FragmentationTerminus.N, 4, 4, 1, 0);
            var decoyPeptideTheorProducts = new List<Product> { a, b, c, d };
            MatchedFragmentIon aa = new MatchedFragmentIon(a, 1, 1, 1);
            MatchedFragmentIon bb = new MatchedFragmentIon(b, 2, 2, 1);
            MatchedFragmentIon cc = new MatchedFragmentIon(c, 3, 3, 1);
            MatchedFragmentIon dd = new MatchedFragmentIon(d, 4, 4, 1);
            var peaks = new List<MatchedFragmentIon> { aa, bb, cc, dd };
            var librarySpectrum = new LibrarySpectrum("library", 0, 0, peaks, 0);
            
            string name = "library/0";
            Assert.AreEqual(name,librarySpectrum.Name);

            string librarySpectrumOverrideToString = "Name: library/0\nMW: 0\nComment: Parent=0 RT=0\nNum peaks: 4\n1\t0.25\t\"b1^1/0ppm\"\n2\t0.5\t\"b2^1/0ppm\"\n3\t0.75\t\"b3^1/0ppm\"\n4\t1\t\"b4^1/0ppm\"";
            Assert.AreEqual(librarySpectrumOverrideToString, librarySpectrum.ToString());

            string spectralAngleOnTheFly = "N/A";
            Assert.AreEqual(spectralAngleOnTheFly,librarySpectrum.CalculateSpectralAngleOnTheFly(peaks));
        }
    }
}
