using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;

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

        [Test]
        public static void CrosslinkPsmFromTsvTest()
        {
            string psmFile = @"FileReadingTests\SearchResults\XL_Intralinks.tsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(6, parsedPsms.Count);
            Assert.That(parsedPsms[0].UniqueSequence, Is.EqualTo("LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)"));
        }

        [Test]
        public static void CrosslinkPsmFromTsvToLibrarySpectrumTest()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\XL_Intralinks.tsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out warnings).ToList();
            Assert.That(warnings.Count == 0);

            CrosslinkLibrarySpectrum librarySpectrum = psms[0].ToLibrarySpectrum() as CrosslinkLibrarySpectrum;
            Assert.IsNotNull(librarySpectrum);
            Assert.AreEqual("Name: LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)/5", librarySpectrum.ToString().Split('\n')[0].Trim());

            // This test would be better if MatchedIon.equals method worked, but it breaks because the mz comparison is implemented incorrectly.
            CollectionAssert.AreEquivalent(librarySpectrum.MatchedFragmentIons.Select(ion => ion.Annotation), psms[0].MatchedIons.Select(ion => ion.Annotation));
            CollectionAssert.AreEquivalent(librarySpectrum.BetaPeptideSpectrum.MatchedFragmentIons.Select(ion => ion.Annotation), psms[0].BetaPeptideMatchedIons.Select(ion => ion.Annotation));
        }
    }
}
