using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework.Legacy;
using Readers;
using MzLibUtil;

namespace Test.FileReadingTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class TestLightweightPsmFromTsv
    {
        [Test]
        [TestCase("oglycoSinglePsms.psmtsv", 2)]
        [TestCase("oGlycoAllPsms.psmtsv", 10)]
        [TestCase("nglyco_f5.psmtsv", 5)]
        [TestCase("nglyco_f5_NewVersion.psmtsv", 5)]
        [TestCase("VariantCrossTest.psmtsv", 15)]
        [TestCase("XL_Intralinks.tsv", 6)]
        [TestCase("XLink.psmtsv", 19)]
        public static void TestLightweightReaderWithMultipleEntryPoints(string path, int expected)
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults", path);

            // Read via lightweight reader
            var lightweightResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out var lwWarnings);
            Assert.AreEqual(expected, lightweightResults.Count);

            // Read via FileReader
            var lwFile = FileReader.ReadFile<LightWeightSpectralMatchFile>(psmFilePath);
            Assert.AreEqual(expected, lwFile.Results.Count);

            // Verify shared fields match between lightweight and full reader
            List<PsmFromTsv> fullResults = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out _);
            Assert.AreEqual(fullResults.Count, lightweightResults.Count);
            for (int i = 0; i < fullResults.Count; i++)
            {
                Assert.AreEqual(fullResults[i].FullSequence, lightweightResults[i].FullSequence);
                Assert.AreEqual(fullResults[i].BaseSeq, lightweightResults[i].BaseSequence);
                Assert.AreEqual(fullResults[i].Ms2ScanNumber, lightweightResults[i].OneBasedScanNumber);
                Assert.AreEqual(fullResults[i].Score, lightweightResults[i].Score);
                Assert.AreEqual(fullResults[i].QValue, lightweightResults[i].QValue);
                Assert.AreEqual(fullResults[i].IsDecoy, lightweightResults[i].IsDecoy);
                Assert.AreEqual(fullResults[i].FileNameWithoutExtension, lightweightResults[i].FileName);
            }
        }

        [Test]
        public static void TestLightweightFileReader()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");

            // Full reader for comparison
            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out var fullWarnings);

            // FileReader.ReadFile<LightWeightSpectralMatchFile>
            var lwFile = FileReader.ReadFile<LightWeightSpectralMatchFile>(psmTsvPath);
            Assert.AreEqual(fullPsms.Count, lwFile.Results.Count);
            Assert.AreEqual(SupportedFileType.psmtsv, lwFile.FileType);
            Assert.Throws<NotImplementedException>(() => lwFile.WriteResults(""));
        }

        [Test]
        public static void TestLightweightDirectConstruction()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");

            var lwFile = new LightWeightSpectralMatchFile(psmTsvPath);
            lwFile.LoadResults();

            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out _);
            Assert.AreEqual(fullPsms.Count, lwFile.Results.Count);
        }

        [Test]
        public static void TestLightweightSharedFieldsMatchFull()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\TDGPTMDSearchResults.psmtsv");

            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out _);
            var lwResults = LightWeightSpectralMatchReader.ReadTsv(psmTsvPath, out _);

            Assert.AreEqual(fullPsms.Count, lwResults.Count);
            for (int i = 0; i < fullPsms.Count; i++)
            {
                var full = fullPsms[i];
                var lw = lwResults[i];

                // ISpectralMatch fields
                Assert.AreEqual(full.FullSequence, lw.FullSequence);
                Assert.AreEqual(full.BaseSeq, lw.BaseSequence);
                Assert.AreEqual(full.Ms2ScanNumber, lw.OneBasedScanNumber);
                Assert.AreEqual(full.Score, lw.Score);
                Assert.AreEqual(full.IsDecoy, lw.IsDecoy);
                Assert.AreEqual(full.FileNameWithoutExtension, lw.FileName);
                Assert.AreEqual(full.FileNameWithoutExtension, lw.FullFilePath);

                // IQuantifiableRecord fields
                Assert.AreEqual(full.PrecursorCharge, lw.ChargeState);
                Assert.AreEqual(full.RetentionTime, lw.RetentionTime);
                Assert.AreEqual(full.MonoisotopicMass, lw.MonoisotopicMass);
                Assert.AreEqual(full.QValue, lw.QValue);
                Assert.AreEqual(full.PEP_QValue, lw.PepQValue);

                // Accession
                Assert.AreEqual(full.Accession, lw.Accession);
            }
        }

        [Test]
        public static void TestLightweightProteinGroupInfos()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out _);
            var lwResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);
            Assert.AreEqual(fullPsms.Count, lwResults.Count);

            // Find an ambiguous result to test ProteinGroupInfos
            var fullAmbiguous = fullPsms.First(p => p.Accession.Contains("|"));
            var lwAmbiguous = lwResults.First(p => p.Accession.Contains("|"));

            Assert.AreEqual(fullAmbiguous.Accession, lwAmbiguous.Accession);
            Assert.AreEqual(fullAmbiguous.ProteinGroupInfos.Count, lwAmbiguous.ProteinGroupInfos.Count);

            for (int j = 0; j < fullAmbiguous.ProteinGroupInfos.Count; j++)
            {
                Assert.AreEqual(fullAmbiguous.ProteinGroupInfos[j].proteinAccessions,
                    lwAmbiguous.ProteinGroupInfos[j].proteinAccessions);
                Assert.AreEqual(fullAmbiguous.ProteinGroupInfos[j].geneName,
                    lwAmbiguous.ProteinGroupInfos[j].geneName);
                Assert.AreEqual(fullAmbiguous.ProteinGroupInfos[j].organism,
                    lwAmbiguous.ProteinGroupInfos[j].organism);
            }
        }

        [Test]
        [TestCase("BottomUpExample.psmtsv")]
        [TestCase("TDGPTMDSearchResults.psmtsv")]
        public static void TestLightweightProteinGroupInfosMultipleFiles(string path)
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults", path);

            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFilePath, out _);
            var lwResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);

            var fullAmbiguous = fullPsms.First(p => p.Accession.Contains("|"));
            var lwAmbiguous = lwResults.First(p => p.Accession.Contains("|"));
            int ambiguityCount = fullAmbiguous.Accession.Count(c => c == '|') + 1;
            Assert.AreEqual(ambiguityCount, lwAmbiguous.ProteinGroupInfos.Count);

            CollectionAssert.AreEquivalent(
                fullAmbiguous.ProteinGroupInfos.Select(g => g.proteinAccessions).ToList(),
                lwAmbiguous.ProteinGroupInfos.Select(g => g.proteinAccessions).ToList());
        }

        [Test]
        public static void TestLightweightGetIdentifiedBioPolymersThrows()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            var lwResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);
            Assert.That(lwResults.Count > 0);
            Assert.Throws<NotSupportedException>(() =>
            {
                lwResults[0].GetIdentifiedBioPolymersWithSetMods();
            });
        }

        [Test]
        public static void TestLightweightCompareTo()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            var lwResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);
            Assert.That(lwResults.Count >= 2);

            // Sort by score descending (CompareTo uses higher-is-better convention)
            var sorted = lwResults.OrderBy(r => r).ToList();
            for (int i = 0; i < sorted.Count - 1; i++)
            {
                Assert.That(sorted[i].Score >= sorted[i + 1].Score);
            }
        }

        [Test]
        public static void TestLightweightRowFilter()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            // Row filter: only targets
            var rowFilters = new Dictionary<string, Func<string, bool>>
            {
                { SpectrumMatchFromTsvHeader.DecoyContaminantTarget, val => val.Trim() == "T" }
            };

            var allResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);
            var filteredResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _, rowFilters: rowFilters);

            Assert.That(filteredResults.Count > 0);
            Assert.That(filteredResults.Count <= allResults.Count);
            Assert.That(filteredResults.All(r => !r.IsDecoy));
        }

        [Test]
        public static void TestLightweightTerminatingFilter()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            double qValueCutoff = 0.01;
            var terminatingFilters = new Dictionary<string, Func<string, bool>>
            {
                { SpectrumMatchFromTsvHeader.QValue, val =>
                    double.TryParse(val.Trim(), NumberStyles.Any, CultureInfo.InvariantCulture, out double q) && q <= qValueCutoff }
            };

            var allResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);
            var filteredResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _,
                terminatingFilters: terminatingFilters);

            Assert.That(filteredResults.Count > 0);
            Assert.That(filteredResults.Count <= allResults.Count);
            Assert.That(filteredResults.All(r => r.QValue <= qValueCutoff));
        }

        [Test]
        public static void TestLightweightCombinedFilters()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            double qValueCutoff = 0.01;
            var rowFilters = new Dictionary<string, Func<string, bool>>
            {
                { SpectrumMatchFromTsvHeader.DecoyContaminantTarget, val => val.Trim() == "T" }
            };
            var terminatingFilters = new Dictionary<string, Func<string, bool>>
            {
                { SpectrumMatchFromTsvHeader.QValue, val =>
                    double.TryParse(val.Trim(), NumberStyles.Any, CultureInfo.InvariantCulture, out double q) && q <= qValueCutoff }
            };

            var filteredResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _,
                rowFilters: rowFilters, terminatingFilters: terminatingFilters);

            Assert.That(filteredResults.Count > 0);
            Assert.That(filteredResults.All(r => !r.IsDecoy));
            Assert.That(filteredResults.All(r => r.QValue <= qValueCutoff));
        }

        [Test]
        public static void TestLightweightWithNewHeaderTerms()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\NewHeader.psmtsv");

            var lwFile = FileReader.ReadFile<LightWeightSpectralMatchFile>(psmTsvPath);
            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out _);
            Assert.AreEqual(fullPsms.Count, lwFile.Results.Count);

            var testResult = lwFile.Results.First();
            Assert.That(!double.IsNaN(testResult.PepQValue));
            Assert.That(!double.IsNaN(testResult.RetentionTime) && testResult.RetentionTime != -1);
        }

        [Test]
        public static void TestLightweightXLink()
        {
            string psmTsvPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\XLink.psmtsv");

            var lwFile = FileReader.ReadFile<LightWeightSpectralMatchFile>(psmTsvPath);
            List<PsmFromTsv> fullPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmTsvPath, out _);
            Assert.AreEqual(fullPsms.Count, lwFile.Results.Count);
        }

        [Test]
        [TestCase("BottomUpExample.psmtsv", "04-30-13_CAST_Frac5_4uL.raw", @"D:/Data/This/Is/A/Folder/Tree/04-30-13_CAST_Frac4_6uL.mzML")]
        [TestCase("TDGPTMDSearchResults.psmtsv", "TDGPTMDSearchSingleSpectra.raw", @"C:/Data/FXN3_tr1_032017-calib.mzML")]
        public static void TestLightweightFileNameFilePathDictionary(string path, string filePathA, string filePathB)
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults", path);

            var lwFile = FileReader.ReadFile<LightWeightSpectralMatchFile>(psmFilePath);
            List<string> filePaths = new List<string> { filePathA, filePathB };

            var dictionary = lwFile.FileNameToFilePath(filePaths);
            Assert.That(dictionary.Count == 2);
            Assert.That(lwFile.GetQuantifiableResults().Any(p => p.FileName.Equals(dictionary.Keys.First())));
            Assert.That(lwFile.GetQuantifiableResults().Any(p => p.FileName.Equals(dictionary.Keys.Last())));
        }

        [Test]
        public static void TestLightweightQuantifiableResults()
        {
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            var lwFile = FileReader.ReadFile<LightWeightSpectralMatchFile>(psmFilePath);
            var quantResults = lwFile.GetQuantifiableResults().ToList();
            Assert.AreEqual(lwFile.Results.Count, quantResults.Count);

            // Verify IQuantifiableRecord properties are populated
            foreach (var result in quantResults)
            {
                Assert.That(!string.IsNullOrEmpty(result.FileName));
                Assert.That(!string.IsNullOrEmpty(result.BaseSequence));
                Assert.That(!string.IsNullOrEmpty(result.FullSequence));
                Assert.That(result.ChargeState > 0);
            }
        }

        [Test]
        public static void TestLightweightIntensitiesNullWhenNoTmt()
        {
            // Standard psmtsv files without TMT columns should have null Intensities
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            var lwResults = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out _);
            Assert.That(lwResults.Count > 0);
            Assert.That(lwResults.All(r => r.Intensities == null));
        }

        [Test]
        public static void TestLightweightFilterOnMissingColumn()
        {
            // Filter on a column that doesn't exist in the file — should not throw
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\SearchResults\BottomUpExample.psmtsv");

            var rowFilters = new Dictionary<string, Func<string, bool>>
            {
                { "NonExistentColumn", _ => true }
            };

            var results = LightWeightSpectralMatchReader.ReadTsv(psmFilePath, out var warnings, rowFilters: rowFilters);
            Assert.That(results.Count > 0);
        }
    }
}
