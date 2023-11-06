using NUnit.Framework;
using Proteomics.Fragmentation;
using Proteomics.PSM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
            MatchedFragmentIon aa = new MatchedFragmentIon(ref a, 1, 1, 1);
            MatchedFragmentIon bb = new MatchedFragmentIon(ref b, 2, 2, 1);
            MatchedFragmentIon cc = new MatchedFragmentIon(ref c, 3, 3, 1);
            MatchedFragmentIon dd = new MatchedFragmentIon(ref d, 4, 4, 1);
            var peaks = new List<MatchedFragmentIon> { aa, bb, cc, dd };
            var librarySpectrum = new LibrarySpectrum("library", 0, 0, peaks, 0);
            //Assert.That(librarySpectrum. .NeutralTheoreticalProduct.ProductType == ProductType.b && decoySpectum[0].NeutralTheoreticalProduct.FragmentNumber == 1 && decoySpectum[0].Intensity == 1);
            //Assert.That(decoySpectum[1].NeutralTheoreticalProduct.ProductType == ProductType.b && decoySpectum[1].NeutralTheoreticalProduct.FragmentNumber == 2 && decoySpectum[1].Intensity == 2);
            //Assert.That(decoySpectum[2].NeutralTheoreticalProduct.ProductType == ProductType.b && decoySpectum[2].NeutralTheoreticalProduct.FragmentNumber == 3 && decoySpectum[2].Intensity == 3);
            //Assert.That(decoySpectum[3].NeutralTheoreticalProduct.ProductType == ProductType.b && decoySpectum[3].NeutralTheoreticalProduct.FragmentNumber == 4 && decoySpectum[3].Intensity == 4);
        }
    }
}
