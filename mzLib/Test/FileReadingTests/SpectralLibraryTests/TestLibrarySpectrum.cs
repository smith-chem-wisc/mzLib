using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using Omics.SpectralLibrary;
using Omics.SpectrumMatch;
using Readers;
using Readers.SpectralLibrary;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.FileReadingTests.SpectralLibraryTests
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
            NUnit.Framework.Assert.That(mspString, Does.Contain("bIb[3-6]"),
                "MSP output should contain internal fragment annotation bIb[3-6]");
            NUnit.Framework.Assert.That(mspString, Does.Contain("bIb[2-5]"),
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

            NUnit.Framework.Assert.That(ex.Message, Does.Contain("Unable to parse fragment ion annotation"));
            NUnit.Framework.Assert.That(ex.Message, Does.Contain("xyz123abc"));
        }

        [Test]
        public static void CrosslinkPsmFromTsvTest()
        {
            string psmFile = @"FileReadingTests\SearchResults\XL_Intralinks.tsv";
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            Assert.AreEqual(6, parsedPsms.Count);
            NUnit.Framework.Assert.That(parsedPsms[0].UniqueSequence, Is.EqualTo("LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)"));
        }

        [Test]
        public static void CrosslinkPsmFromTsvToLibrarySpectrumTest()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults\XL_Intralinks.tsv");
            List<string> warnings = new();
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out warnings).ToList();
            NUnit.Framework.Assert.That(warnings.Count == 0);

            CrosslinkLibrarySpectrum librarySpectrum = psms[0].ToLibrarySpectrum() as CrosslinkLibrarySpectrum;
            Assert.IsNotNull(librarySpectrum);
            Assert.AreEqual("Name: LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)/5", librarySpectrum.ToString().Split('\n')[0].Trim());

            // This test would be better if MatchedIon.equals method worked, but it breaks because the mz comparison is implemented incorrectly.
            CollectionAssert.AreEquivalent(librarySpectrum.MatchedFragmentIons.Select(ion => ion.Annotation), psms[0].MatchedIons.Select(ion => ion.Annotation));
            CollectionAssert.AreEquivalent(librarySpectrum.BetaPeptideSpectrum.MatchedFragmentIons.Select(ion => ion.Annotation), psms[0].BetaPeptideMatchedIons.Select(ion => ion.Annotation));
        }

        [Test]
        public static void LibrarySpectrum_Equals_Null_ReturnsFalse()
        {
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 1000, 1) };
            var spectrum = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            Assert.IsFalse(spectrum.Equals((LibrarySpectrum?)null));
        }

        [Test]
        public static void LibrarySpectrum_Equals_SameReference_ReturnsTrue()
        {
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 1000, 1) };
            var spectrum = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            Assert.IsTrue(spectrum.Equals(spectrum));
        }

        [Test]
        public static void LibrarySpectrum_Equals_SameValues_ReturnsTrue()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var peaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks1, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks2, 10.0);
            Assert.IsTrue(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_DifferentSequence_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            var spectrum2 = new LibrarySpectrum("PEP", 500, 2, peaks, 10.0);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_DifferentRetentionTime_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks, 20.0);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_DifferentPrecursorMz_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 501, 2, peaks, 10.0);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_DifferentChargeState_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 3, peaks, 10.0);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_DifferentMatchedFragmentIons_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var peaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 101, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks1, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks2, 10.0);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_SameArraysDifferentAnnotations_ReturnsFalse()
        {
            var product1 = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var product2 = new Product(ProductType.y, FragmentationTerminus.C, 1, 2, 2, 0);
            var peaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product1, 100, 1000, 1) };
            var peaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product2, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks1, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks2, 10.0);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_DifferentIsDecoy_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0, isDecoy: false);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0, isDecoy: true);
            Assert.IsFalse(spectrum1.Equals(spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_ObjectEquals_Null_ReturnsFalse()
        {
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 1000, 1) };
            var spectrum = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            Assert.IsFalse(spectrum.Equals((object?)null));
        }

        [Test]
        public static void LibrarySpectrum_ObjectEquals_WrongType_ReturnsFalse()
        {
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 100, 1000, 1) };
            var spectrum = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            Assert.IsFalse(spectrum.Equals("not a spectrum"));
        }

        [Test]
        public static void LibrarySpectrum_ObjectEquals_SameValues_ReturnsTrue()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var peaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks1, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks2, 10.0);
            Assert.IsTrue(spectrum1.Equals((object)spectrum2));
        }

        [Test]
        public static void LibrarySpectrum_Equals_CrosslinkLibrarySpectrum_ReturnsFalse()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var librarySpectrum = new LibrarySpectrum("SEQ", 500, 2, peaks, 10.0);
            var crosslinkSpectrum = new CrosslinkLibrarySpectrum("SEQ(1)PEP", 500, 2, peaks, 10.0, peaks);

            Assert.IsFalse(librarySpectrum.Equals(crosslinkSpectrum));
            Assert.IsFalse(librarySpectrum.Equals((object)crosslinkSpectrum));
        }

        [Test]
        public static void CrosslinkLibrarySpectrum_Equals_SameValues_ReturnsTrue()
        {
            var alphaProduct = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var betaProduct = new Product(ProductType.y, FragmentationTerminus.C, 1, 1, 1, 0);

            var alphaPeaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(alphaProduct, 100, 1000, 1) };
            var betaPeaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(betaProduct, 200, 500, 1) };
            var alphaPeaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(alphaProduct, 100, 1000, 1) };
            var betaPeaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(betaProduct, 200, 500, 1) };

            var spectrum1 = new CrosslinkLibrarySpectrum("SEQ(1)PEP", 500, 2, alphaPeaks1, 10.0, betaPeaks1);
            var spectrum2 = new CrosslinkLibrarySpectrum("SEQ(1)PEP", 500, 2, alphaPeaks2, 10.0, betaPeaks2);

            Assert.IsTrue(spectrum1.Equals((object)spectrum2));
            Assert.AreEqual(spectrum1.GetHashCode(), spectrum2.GetHashCode());
        }

        [Test]
        public static void LibrarySpectrum_GetHashCode_EqualObjects_SameHashCode()
        {
            var product = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            var peaks1 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var peaks2 = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, 100, 1000, 1) };
            var spectrum1 = new LibrarySpectrum("SEQ", 500, 2, peaks1, 10.0);
            var spectrum2 = new LibrarySpectrum("SEQ", 500, 2, peaks2, 10.0);
            Assert.AreEqual(spectrum1.GetHashCode(), spectrum2.GetHashCode());
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToLibrary_IdenticalSpectra()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib1 = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var lib2 = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);

            var similarity = lib1.ComputeSpectralSimilarity(lib2, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToLibrary_DifferentSpectra()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks1 = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var peaks2 = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 600.0, 3000.0, 1),
                new MatchedFragmentIon(b, 700.0, 4000.0, 1)
            };
            var lib1 = new LibrarySpectrum("SEQ", 500.0, 2, peaks1, 10.0);
            var lib2 = new LibrarySpectrum("SEQ", 500.0, 2, peaks2, 10.0);

            var similarity = lib1.ComputeSpectralSimilarity(lib2, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(0.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToMzSpectrum()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var mz = new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);

            var similarity = lib.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_MzSpectrumToLibrary()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var mz = new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);

            var similarity = mz.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_PartialOverlap()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks1 = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var peaks2 = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 600.0, 3000.0, 1)
            };
            var lib1 = new LibrarySpectrum("SEQ", 500.0, 2, peaks1, 10.0);
            var lib2 = new LibrarySpectrum("SEQ", 500.0, 2, peaks2, 10.0);

            var similarity = lib1.ComputeSpectralSimilarity(lib2, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            // Partial overlap should yield a cosine between 0 and 1
            double? cosine = similarity.CosineSimilarity();
            NUnit.Framework.Assert.That(cosine, Is.GreaterThan(0.0).And.LessThan(1.0));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToArrays()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            double[] xArray = new double[] { 400.0, 500.0 };
            double[] yArray = new double[] { 1000.0, 2000.0 };

            var similarity = lib.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_ArraysToLibrary()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            double[] xArray = new double[] { 400.0, 500.0 };
            double[] yArray = new double[] { 1000.0, 2000.0 };

            var similarity = xArray.ComputeSpectralSimilarity(yArray, lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        #region MsDataScan overload tests

        private static MsDataScan CreateTestScan(double[] x, double[] y)
        {
            var spectrum = new MzSpectrum(x, y, false);
            return new MsDataScan(
                spectrum, 1, 1, true, Polarity.Positive, 1.0,
                new MzRange(300, 2000), "f", MZAnalyzerType.Orbitrap,
                spectrum.SumOfAllY, null, null, "scan=1");
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_MsDataScanToLibrary()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var similarity = scan.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_MsDataScanToMzSpectrum()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var mz = new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);

            var similarity = scan.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToMsDataScan()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var similarity = lib.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_MzSpectrumToMsDataScan()
        {
            var mz = new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var similarity = mz.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_MsDataScanToArrays()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            double[] xArray = new double[] { 400.0, 500.0 };
            double[] yArray = new double[] { 1000.0, 2000.0 };

            var similarity = scan.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_ArraysToMsDataScan()
        {
            double[] xArray = new double[] { 400.0, 500.0 };
            double[] yArray = new double[] { 1000.0, 2000.0 };
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var similarity = xArray.ComputeSpectralSimilarity(yArray, scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.IsNotNull(similarity);
            NUnit.Framework.Assert.That(similarity.CosineSimilarity(), Is.EqualTo(1.0).Within(0.001));
        }

        #endregion

        #region SpectralSimilarity dispatcher tests

        [Test]
        public static void SpectralSimilarity_GetSimilarityMeasure_MatchesDirectCalls()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var similarity = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            // Verify each enum-backed call returns the same value as the direct method call
            Assert.AreEqual(similarity.CosineSimilarity(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.CosineSimilarity), "CosineSimilarity mismatch");
            Assert.AreEqual(similarity.SpectralContrastAngle(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.SpectralContrastAngle), "SpectralContrastAngle mismatch");
            Assert.AreEqual(similarity.EuclideanDistance(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.EuclideanDistance), "EuclideanDistance mismatch");
            Assert.AreEqual(similarity.BrayCurtis(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.BrayCurtis), "BrayCurtis mismatch");
            Assert.AreEqual(similarity.PearsonsCorrelation(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.PearsonsCorrelation), "PearsonsCorrelation mismatch");
            Assert.AreEqual(similarity.DotProduct(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.DotProduct), "DotProduct mismatch");
            Assert.AreEqual(similarity.SpectralEntropy(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.SpectralEntropy), "SpectralEntropy mismatch");
            Assert.AreEqual(similarity.KullbackLeiblerDivergence_P_Q(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.KullbackLeiblerDivergence_P_Q), "KullbackLeiblerDivergence_P_Q mismatch");
            Assert.AreEqual(similarity.SearleSimilarity(), similarity.GetSimilarityMeasure(SpectralSimilarity.SimilarityMeasures.SearleSimilarity), "SearleSimilarity mismatch");
        }

        [Test]
        public static void SpectralSimilarity_GetAllSimilarityMeasures_ReturnsAllMeasures()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var similarity = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var all = similarity.GetAllSimilarityMeasures().ToList();
            int expectedCount = Enum.GetValues(typeof(SpectralSimilarity.SimilarityMeasures)).Length;
            Assert.AreEqual(expectedCount, all.Count);

            foreach (var (measure, value) in all)
            {
                Assert.AreEqual(similarity.GetSimilarityMeasure(measure), value, $"Mismatch for {measure}");
            }
        }

        [Test]
        public static void SpectralSimilarity_GetAllSimilarityMeasuresExcept_ExcludesSpecified()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var similarity = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var except = similarity.GetAllSimilarityMeasuresExcept(SpectralSimilarity.SimilarityMeasures.CosineSimilarity).ToList();
            int expectedCount = Enum.GetValues(typeof(SpectralSimilarity.SimilarityMeasures)).Length - 1;
            Assert.AreEqual(expectedCount, except.Count);
            Assert.IsFalse(except.Any(m => m.Item1 == SpectralSimilarity.SimilarityMeasures.CosineSimilarity));
        }

        [Test]
        public static void SpectralSimilarity_GetSelectedSimilarityMeasures_ReturnsOnlySelected()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var similarity = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var selected = similarity.GetSelectedSimilarityMeasures(
                SpectralSimilarity.SimilarityMeasures.CosineSimilarity,
                SpectralSimilarity.SimilarityMeasures.DotProduct).ToList();

            Assert.AreEqual(2, selected.Count);
            Assert.IsTrue(selected.Any(m => m.Item1 == SpectralSimilarity.SimilarityMeasures.CosineSimilarity));
            Assert.IsTrue(selected.Any(m => m.Item1 == SpectralSimilarity.SimilarityMeasures.DotProduct));
        }

        [Test]
        public static void SpectralSimilarity_GetSimilarityMeasure_InvalidEnum_ReturnsNull()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var similarity = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var invalid = (SpectralSimilarity.SimilarityMeasures)999;
            Assert.IsNull(similarity.GetSimilarityMeasure(invalid));
        }

        #endregion

        #region Extension shortcut tests

        [Test]
        public static void LibrarySpectrum_GetSimilarityMeasure_ShortcutMatchesExplicit()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);

            var shortcut = lib.GetSimilarityMeasure(lib, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetSimilarityMeasure_MzSpectrum_ShortcutMatchesExplicit()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var mz = new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);

            var shortcut = lib.GetSimilarityMeasure(mz, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = lib.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MzSpectrum_GetSimilarityMeasure_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var mz = new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);

            var shortcut = mz.GetSimilarityMeasure(lib, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = mz.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetSimilarityMeasure_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = scan.GetSimilarityMeasure(lib, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = scan.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetAllSimilarityMeasures_ShortcutMatchesExplicit()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);

            var shortcut = lib.GetAllSimilarityMeasures(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            Assert.AreEqual(explicitValue.Count, shortcut.Count);
            for (int i = 0; i < explicitValue.Count; i++)
            {
                Assert.AreEqual(explicitValue[i].Item1, shortcut[i].Item1);
                Assert.AreEqual(explicitValue[i].Item2, shortcut[i].Item2);
            }
        }

        [Test]
        public static void LibrarySpectrum_GetSelectedSimilarityMeasures_ShortcutMatchesExplicit()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            var peaks = new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, peaks, 10.0);
            var measures = new[] { SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SimilarityMeasures.DotProduct };

            var shortcut = lib.GetSelectedSimilarityMeasures(lib, measures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(measures).ToList();

            Assert.AreEqual(explicitValue.Count, shortcut.Count);
            for (int i = 0; i < explicitValue.Count; i++)
            {
                Assert.AreEqual(explicitValue[i].Item1, shortcut[i].Item1);
                Assert.AreEqual(explicitValue[i].Item2, shortcut[i].Item2);
            }
        }

        #endregion

        #region Shared fixture helpers

        private static List<MatchedFragmentIon> CreateTestPeaks()
        {
            Product a = new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0);
            Product b = new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0);
            return new List<MatchedFragmentIon>
            {
                new MatchedFragmentIon(a, 400.0, 1000.0, 1),
                new MatchedFragmentIon(b, 500.0, 2000.0, 1)
            };
        }

        private static LibrarySpectrum CreateTestLibrarySpectrum()
        {
            return new LibrarySpectrum("SEQ", 500.0, 2, CreateTestPeaks(), 10.0);
        }

        private static MzSpectrum CreateTestMzSpectrum()
        {
            return new MzSpectrum(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 }, false);
        }

        private static (double[] xArray, double[] yArray) CreateTestArrays()
        {
            return (new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
        }

        private static void AssertSimilarityMeasuresEqual(
            List<(SpectralSimilarity.SimilarityMeasures, double?)> expected,
            List<(SpectralSimilarity.SimilarityMeasures, double?)> actual)
        {
            Assert.AreEqual(expected.Count, actual.Count);
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.AreEqual(expected[i].Item1, actual[i].Item1);
                Assert.AreEqual(expected[i].Item2, actual[i].Item2);
            }
        }

        #endregion

        #region GetSimilarityMeasure shortcut tests — remaining shapes

        [Test]
        public static void LibrarySpectrum_GetSimilarityMeasure_Arrays_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var (xArray, yArray) = CreateTestArrays();

            var shortcut = lib.GetSimilarityMeasure(xArray, yArray, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = lib.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void Arrays_GetSimilarityMeasure_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var (xArray, yArray) = CreateTestArrays();
            var lib = CreateTestLibrarySpectrum();

            var shortcut = xArray.GetSimilarityMeasure(yArray, lib, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = xArray.ComputeSpectralSimilarity(yArray, lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetSimilarityMeasure_MzSpectrum_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var mz = CreateTestMzSpectrum();

            var shortcut = scan.GetSimilarityMeasure(mz, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = scan.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetSimilarityMeasure_MsDataScan_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = lib.GetSimilarityMeasure(scan, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = lib.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MzSpectrum_GetSimilarityMeasure_MsDataScan_ShortcutMatchesExplicit()
        {
            var mz = CreateTestMzSpectrum();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = mz.GetSimilarityMeasure(scan, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = mz.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetSimilarityMeasure_Arrays_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var (xArray, yArray) = CreateTestArrays();

            var shortcut = scan.GetSimilarityMeasure(xArray, yArray, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = scan.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        [Test]
        public static void Arrays_GetSimilarityMeasure_MsDataScan_ShortcutMatchesExplicit()
        {
            var (xArray, yArray) = CreateTestArrays();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = xArray.GetSimilarityMeasure(yArray, scan, SpectralSimilarity.SimilarityMeasures.CosineSimilarity, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);
            var explicitValue = xArray.ComputeSpectralSimilarity(yArray, scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).CosineSimilarity();

            Assert.AreEqual(explicitValue, shortcut);
        }

        #endregion

        #region GetAllSimilarityMeasures shortcut tests — remaining shapes

        [Test]
        public static void LibrarySpectrum_GetAllSimilarityMeasures_MzSpectrum_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var mz = CreateTestMzSpectrum();

            var shortcut = lib.GetAllSimilarityMeasures(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MzSpectrum_GetAllSimilarityMeasures_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var mz = CreateTestMzSpectrum();
            var lib = CreateTestLibrarySpectrum();

            var shortcut = mz.GetAllSimilarityMeasures(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = mz.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetAllSimilarityMeasures_Arrays_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var (xArray, yArray) = CreateTestArrays();

            var shortcut = lib.GetAllSimilarityMeasures(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void Arrays_GetAllSimilarityMeasures_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var (xArray, yArray) = CreateTestArrays();
            var lib = CreateTestLibrarySpectrum();

            var shortcut = xArray.GetAllSimilarityMeasures(yArray, lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = xArray.ComputeSpectralSimilarity(yArray, lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetAllSimilarityMeasures_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var lib = CreateTestLibrarySpectrum();

            var shortcut = scan.GetAllSimilarityMeasures(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = scan.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetAllSimilarityMeasures_MzSpectrum_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var mz = CreateTestMzSpectrum();

            var shortcut = scan.GetAllSimilarityMeasures(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = scan.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetAllSimilarityMeasures_MsDataScan_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = lib.GetAllSimilarityMeasures(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MzSpectrum_GetAllSimilarityMeasures_MsDataScan_ShortcutMatchesExplicit()
        {
            var mz = CreateTestMzSpectrum();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = mz.GetAllSimilarityMeasures(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = mz.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetAllSimilarityMeasures_Arrays_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var (xArray, yArray) = CreateTestArrays();

            var shortcut = scan.GetAllSimilarityMeasures(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = scan.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void Arrays_GetAllSimilarityMeasures_MsDataScan_ShortcutMatchesExplicit()
        {
            var (xArray, yArray) = CreateTestArrays();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = xArray.GetAllSimilarityMeasures(yArray, scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = xArray.ComputeSpectralSimilarity(yArray, scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetAllSimilarityMeasures().ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        #endregion

        #region GetSelectedSimilarityMeasures shortcut tests — remaining shapes

        private static readonly SpectralSimilarity.SimilarityMeasures[] _selectedMeasures =
        {
            SpectralSimilarity.SimilarityMeasures.CosineSimilarity,
            SpectralSimilarity.SimilarityMeasures.DotProduct
        };

        [Test]
        public static void LibrarySpectrum_GetSelectedSimilarityMeasures_MzSpectrum_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var mz = CreateTestMzSpectrum();

            var shortcut = lib.GetSelectedSimilarityMeasures(mz, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MzSpectrum_GetSelectedSimilarityMeasures_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var mz = CreateTestMzSpectrum();
            var lib = CreateTestLibrarySpectrum();

            var shortcut = mz.GetSelectedSimilarityMeasures(lib, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = mz.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetSelectedSimilarityMeasures_Arrays_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var (xArray, yArray) = CreateTestArrays();

            var shortcut = lib.GetSelectedSimilarityMeasures(xArray, yArray, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void Arrays_GetSelectedSimilarityMeasures_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var (xArray, yArray) = CreateTestArrays();
            var lib = CreateTestLibrarySpectrum();

            var shortcut = xArray.GetSelectedSimilarityMeasures(yArray, lib, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = xArray.ComputeSpectralSimilarity(yArray, lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetSelectedSimilarityMeasures_LibrarySpectrum_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var lib = CreateTestLibrarySpectrum();

            var shortcut = scan.GetSelectedSimilarityMeasures(lib, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = scan.ComputeSpectralSimilarity(lib, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetSelectedSimilarityMeasures_MzSpectrum_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var mz = CreateTestMzSpectrum();

            var shortcut = scan.GetSelectedSimilarityMeasures(mz, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = scan.ComputeSpectralSimilarity(mz, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void LibrarySpectrum_GetSelectedSimilarityMeasures_MsDataScan_ShortcutMatchesExplicit()
        {
            var lib = CreateTestLibrarySpectrum();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = lib.GetSelectedSimilarityMeasures(scan, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = lib.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MzSpectrum_GetSelectedSimilarityMeasures_MsDataScan_ShortcutMatchesExplicit()
        {
            var mz = CreateTestMzSpectrum();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = mz.GetSelectedSimilarityMeasures(scan, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = mz.ComputeSpectralSimilarity(scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void MsDataScan_GetSelectedSimilarityMeasures_Arrays_ShortcutMatchesExplicit()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var (xArray, yArray) = CreateTestArrays();

            var shortcut = scan.GetSelectedSimilarityMeasures(xArray, yArray, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = scan.ComputeSpectralSimilarity(xArray, yArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        [Test]
        public static void Arrays_GetSelectedSimilarityMeasures_MsDataScan_ShortcutMatchesExplicit()
        {
            var (xArray, yArray) = CreateTestArrays();
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });

            var shortcut = xArray.GetSelectedSimilarityMeasures(yArray, scan, _selectedMeasures, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).ToList();
            var explicitValue = xArray.ComputeSpectralSimilarity(yArray, scan, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true).GetSelectedSimilarityMeasures(_selectedMeasures).ToList();

            AssertSimilarityMeasuresEqual(explicitValue, shortcut);
        }

        #endregion

        #region Asymmetric similarity tests (catches swapped-argument bugs)

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToMzSpectrum_Asymmetric()
        {
            var expPeaks = new List<MatchedFragmentIon>
            {
                new(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 400.0, 1000.0, 1),
                new(new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0), 500.0, 2000.0, 1),
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, expPeaks, 10.0);
            var theoPeaks = new List<MatchedFragmentIon>
            {
                new(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 400.0, 1000.0, 1),
            };
            var libTheo = new LibrarySpectrum("SEQ", 500.0, 2, theoPeaks, 10.0);

            var result = lib.ComputeSpectralSimilarity((MzSpectrum)libTheo,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var explicitResult = new SpectralSimilarity(
                lib.XArray, lib.YArray,
                libTheo.XArray, libTheo.YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.That(result.CosineSimilarity(), Is.EqualTo(explicitResult.CosineSimilarity()).Within(0.001));
        }

        [Test]
        public static void LibrarySpectrum_ComputeSpectralSimilarity_LibraryToArrays_Asymmetric()
        {
            var expPeaks = new List<MatchedFragmentIon>
            {
                new(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 400.0, 1000.0, 1),
                new(new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0), 500.0, 2000.0, 1),
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, expPeaks, 10.0);
            var theoX = new double[] { 400.0 };
            var theoY = new double[] { 1000.0 };

            var result = lib.ComputeSpectralSimilarity(theoX, theoY,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var explicitResult = new SpectralSimilarity(
                lib.XArray, lib.YArray,
                theoX, theoY,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.That(result.CosineSimilarity(), Is.EqualTo(explicitResult.CosineSimilarity()).Within(0.001));
        }

        [Test]
        public static void MzSpectrum_ComputeSpectralSimilarity_LibrarySpectrum_Asymmetric()
        {
            var expPeaks = new List<MatchedFragmentIon>
            {
                new(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 400.0, 1000.0, 1),
                new(new Product(ProductType.b, FragmentationTerminus.N, 2, 2, 1, 0), 500.0, 2000.0, 1),
            };
            var expLib = new LibrarySpectrum("SEQ", 500.0, 2, expPeaks, 10.0);
            var theoMz = new MzSpectrum(new double[] { 400.0 }, new double[] { 1000.0 }, false);

            var result = expLib.ComputeSpectralSimilarity(theoMz,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var explicitResult = new SpectralSimilarity(
                expLib.XArray, expLib.YArray,
                theoMz.XArray, theoMz.YArray,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.That(result.CosineSimilarity(), Is.EqualTo(explicitResult.CosineSimilarity()).Within(0.001));
        }

        [Test]
        public static void MsDataScan_ComputeSpectralSimilarity_LibrarySpectrum_Asymmetric()
        {
            var scan = CreateTestScan(new double[] { 400.0, 500.0 }, new double[] { 1000.0, 2000.0 });
            var theoPeaks = new List<MatchedFragmentIon>
            {
                new(new Product(ProductType.b, FragmentationTerminus.N, 1, 1, 1, 0), 400.0, 1000.0, 1),
            };
            var lib = new LibrarySpectrum("SEQ", 500.0, 2, theoPeaks, 10.0);

            var result = scan.ComputeSpectralSimilarity(lib,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            var explicitResult = new SpectralSimilarity(
                scan.MassSpectrum, lib,
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak, 20, true);

            Assert.That(result.CosineSimilarity(), Is.EqualTo(explicitResult.CosineSimilarity()).Within(0.001));
        }

        #endregion
    }
}
