using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;
using Readers.SpectralLibrary;

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
            Assert.AreEqual(name, librarySpectrum.Name);

            string librarySpectrumOverrideToString = "Name: library/0\nMW: 0\nComment: Parent=0 RT=0\nNum peaks: 4\n1\t0.25\t\"b1^1/0ppm\"\n2\t0.5\t\"b2^1/0ppm\"\n3\t0.75\t\"b3^1/0ppm\"\n4\t1\t\"b4^1/0ppm\"";
            Assert.AreEqual(librarySpectrumOverrideToString, librarySpectrum.ToString());

            string spectralAngleOnTheFly = "N/A";
            Assert.AreEqual(spectralAngleOnTheFly, librarySpectrum.CalculateSpectralAngleOnTheFly(peaks));
        }

        [Test]
        public static void TestInternalFragmentIonRoundTrip()
        {
            // Create internal fragment products (bIb type - b-ion on both termini)
            // Internal fragment from residue 3 to residue 6
            var internalProduct1 = new Product(
                ProductType.b,
                FragmentationTerminus.None,
                neutralMass: 400.5,
                fragmentNumber: 3,        // start residue
                residuePosition: 3,
                neutralLoss: 0,
                secondaryProductType: ProductType.b,
                secondaryFragmentNumber: 6);  // end residue

            // Internal fragment from residue 2 to residue 5
            var internalProduct2 = new Product(
                ProductType.b,
                FragmentationTerminus.None,
                neutralMass: 350.3,
                fragmentNumber: 2,
                residuePosition: 2,
                neutralLoss: 0,
                secondaryProductType: ProductType.b,
                secondaryFragmentNumber: 5);

            // Create matched fragment ions
            var ion1 = new MatchedFragmentIon(internalProduct1, 401.5, 1.0, 1);
            var ion2 = new MatchedFragmentIon(internalProduct2, 351.3, 0.8, 1);

            // Verify these are recognized as internal fragments
            Assert.IsTrue(ion1.IsInternalFragment, "Ion1 should be an internal fragment");
            Assert.IsTrue(ion2.IsInternalFragment, "Ion2 should be an internal fragment");

            // Create library spectrum
            var peaks = new List<MatchedFragmentIon> { ion1, ion2 };
            var librarySpectrum = new LibrarySpectrum("PEPTIDEK", 500.0, 2, peaks, 100.0);

            // Convert to MSP string format
            string mspString = librarySpectrum.ToString();
            TestContext.WriteLine("MSP Output:\n" + mspString);

            // Verify the annotation format contains internal fragment notation
            Assert.That(mspString, Does.Contain("bIb[3-6]"),
                "MSP output should contain internal fragment annotation bIb[3-6]");
            Assert.That(mspString, Does.Contain("bIb[2-5]"),
                "MSP output should contain internal fragment annotation bIb[2-5]");

            // Now test reading back the fragment ions using ReadFragmentIon
            char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
            char[] neutralLossSplit = new char[] { '-' };

            // Parse the first internal ion line
            string ionLine1 = "401.5\t1\t\"bIb[3-6]^1/0ppm\"";
            var parsedIon1 = SpectralLibrary.ReadFragmentIon(ionLine1, fragmentSplit, neutralLossSplit, "PEPTIDEK");

            Assert.IsTrue(parsedIon1.IsInternalFragment,
                "Parsed ion should be recognized as internal fragment");
            Assert.AreEqual(ProductType.b, parsedIon1.NeutralTheoreticalProduct.ProductType,
                "Primary product type should be b");
            Assert.AreEqual(ProductType.b, parsedIon1.NeutralTheoreticalProduct.SecondaryProductType,
                "Secondary product type should be b");
            Assert.AreEqual(3, parsedIon1.NeutralTheoreticalProduct.FragmentNumber,
                "Start residue should be 3");
            Assert.AreEqual(6, parsedIon1.NeutralTheoreticalProduct.SecondaryFragmentNumber,
                "End residue should be 6");
            Assert.AreEqual(1, parsedIon1.Charge, "Charge should be 1");

            // Parse the second internal ion line
            string ionLine2 = "351.3\t0.8\t\"bIb[2-5]^1/0ppm\"";
            var parsedIon2 = SpectralLibrary.ReadFragmentIon(ionLine2, fragmentSplit, neutralLossSplit, "PEPTIDEK");

            Assert.IsTrue(parsedIon2.IsInternalFragment,
                "Parsed ion should be recognized as internal fragment");
            Assert.AreEqual(2, parsedIon2.NeutralTheoreticalProduct.FragmentNumber,
                "Start residue should be 2");
            Assert.AreEqual(5, parsedIon2.NeutralTheoreticalProduct.SecondaryFragmentNumber,
                "End residue should be 5");
        }

        [Test]
        public static void TestReadFragmentIon_InternalFragment_VariousFormats()
        {
            char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
            char[] neutralLossSplit = new char[] { '-' };

            // Test various internal fragment annotation formats
            var testCases = new[]
            {
                ("401.5\t1\t\"bIb[3-6]^1/0ppm\"", ProductType.b, ProductType.b, 3, 6, 1),
                ("500.0\t0.5\t\"yIy[2-8]^2/0ppm\"", ProductType.y, ProductType.y, 2, 8, 2),
                ("300.0\t0.75\t\"aIa[1-4]^1/0ppm\"", ProductType.a, ProductType.a, 1, 4, 1),
            };

            foreach (var (ionLine, expectedPrimary, expectedSecondary, expectedStart, expectedEnd, expectedCharge) in testCases)
            {
                var parsedIon = SpectralLibrary.ReadFragmentIon(ionLine, fragmentSplit, neutralLossSplit, "TESTPEPTIDE");

                Assert.IsTrue(parsedIon.IsInternalFragment,
                    $"Ion from '{ionLine}' should be internal fragment");
                Assert.AreEqual(expectedPrimary, parsedIon.NeutralTheoreticalProduct.ProductType,
                    $"Primary type mismatch for '{ionLine}'");
                Assert.AreEqual(expectedSecondary, parsedIon.NeutralTheoreticalProduct.SecondaryProductType,
                    $"Secondary type mismatch for '{ionLine}'");
                Assert.AreEqual(expectedStart, parsedIon.NeutralTheoreticalProduct.FragmentNumber,
                    $"Start residue mismatch for '{ionLine}'");
                Assert.AreEqual(expectedEnd, parsedIon.NeutralTheoreticalProduct.SecondaryFragmentNumber,
                    $"End residue mismatch for '{ionLine}'");
                Assert.AreEqual(expectedCharge, parsedIon.Charge,
                    $"Charge mismatch for '{ionLine}'");
            }
        }
        [Test]
        public static void TestReadFragmentIon_InvalidAnnotation_ThrowsInformativeException()
        {
            char[] fragmentSplit = new char[] { '\t', '\"', ')', '/' };
            char[] neutralLossSplit = new char[] { '-' };

            // Test with a malformed annotation that neither regex can parse
            string malformedLine = "401.5\t1\t\"xyz123abc/0ppm\"";

            var ex = Assert.Throws<MzLibUtil.MzLibException>(() =>
                SpectralLibrary.ReadFragmentIon(malformedLine, fragmentSplit, neutralLossSplit, "PEPTIDE"));

            Assert.That(ex.Message, Does.Contain("Unable to parse fragment ion annotation"));
            Assert.That(ex.Message, Does.Contain("xyz123abc"));
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