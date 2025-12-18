using NUnit.Framework;
using Omics;
using Quantification;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.Quantification
{
    [TestFixture]
    public class WriterTests
    {
        /// <summary>
        /// CRITICAL: Verifies that WriteRawData produces a valid TSV file with correct structure.
        /// This is critical because:
        /// 1. File output is the primary deliverable of quantification - if it fails, no results are usable
        /// 2. Header/data column alignment must match exactly or downstream tools will not parse correctly
        /// 3. Both LFQ (single intensity) and TMT (multiple reporters) paths must work correctly
        /// </summary>
        [Test]
        public void WriteRawData_ProducesCorrectTsvStructure_LfqAndTmt()
        {
            // Arrange
            string outputDir = Path.Combine(Path.GetTempPath(), "QuantWriterTest_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(outputDir);

            try
            {
                // LFQ match: single intensity value
                var lfqMatch = new TestSpectralMatch(
                    filePath: "Sample1.raw",
                    fullSequence: "PEPTIDE[+15.995]K",
                    baseSequence: "PEPTIDEK",
                    score: 50,
                    scanNumber: 100)
                {
                    QuantValues = new double[] { 1000000.5 }
                };

                // TMT match: multiple reporter intensities
                var tmtMatch = new TestSpectralMatch(
                    filePath: "Sample2.raw",
                    fullSequence: "ANOTHERPEPTIDE",
                    baseSequence: "ANOTHERPEPTIDE",
                    score: 75,
                    scanNumber: 200)
                {
                    QuantValues = new double[] { 100.0, 200.0, 300.0 }
                };

                var matches = new List<ISpectralMatch> { lfqMatch, tmtMatch };

                // Act
                QuantificationWriter.WriteRawData(matches, outputDir);

                // Assert
                string outputPath = Path.Combine(outputDir, "RawQuantification.tsv");
                Assert.That(File.Exists(outputPath), Is.True, "Output file should be created");

                var lines = File.ReadAllLines(outputPath);
                Assert.That(lines.Length, Is.EqualTo(3), "Should have 1 header + 2 data rows");

                // Verify header has correct columns (TMT case drives column count)
                var headerColumns = lines[0].Split('\t');
                Assert.That(headerColumns[0], Is.EqualTo("FileName"));
                Assert.That(headerColumns[1], Is.EqualTo("Ms2ScanNumber"));
                Assert.That(headerColumns[2], Is.EqualTo("FullSequence"));
                Assert.That(headerColumns[3], Is.EqualTo("BaseSequence"));
                Assert.That(headerColumns.Length, Is.EqualTo(7), "4 base columns + 3 reporter columns");

                // Verify data row structure matches header
                var dataRow1 = lines[1].Split('\t');
                var dataRow2 = lines[2].Split('\t');
                Assert.That(dataRow1.Length, Is.EqualTo(headerColumns.Length), "Data columns must match header");
                Assert.That(dataRow2.Length, Is.EqualTo(headerColumns.Length), "Data columns must match header");

                // Verify key data is written correctly
                Assert.That(lines.Any(l => l.Contains("PEPTIDE[+15.995]K")), Is.True, "Full sequence should be preserved");
                Assert.That(lines.Any(l => l.Contains("100")), Is.True, "Scan number should be written");
            }
            finally
            {
                // Cleanup
                if (Directory.Exists(outputDir))
                    Directory.Delete(outputDir, recursive: true);
            }
        }

        /// <summary>
        /// CRITICAL: Verifies empty/null input doesn't crash and produces no file.
        /// Prevents runtime exceptions in production when no PSMs pass filters.
        /// </summary>
        [Test]
        public void WriteRawData_EmptyInput_DoesNotThrowOrCreateFile()
        {
            string outputDir = Path.Combine(Path.GetTempPath(), "QuantWriterTest_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(outputDir);

            try
            {
                // Act & Assert - should not throw
                Assert.DoesNotThrow(() => QuantificationWriter.WriteRawData(new List<ISpectralMatch>(), outputDir));
                Assert.DoesNotThrow(() => QuantificationWriter.WriteRawData(null, outputDir));

                string outputPath = Path.Combine(outputDir, "RawQuantification.tsv");
                Assert.That(File.Exists(outputPath), Is.False, "No file should be created for empty input");
            }
            finally
            {
                if (Directory.Exists(outputDir))
                    Directory.Delete(outputDir, recursive: true);
            }
        }

        #region Test Helper - Minimal ISpectralMatch implementation

        /// <summary>
        /// Minimal test implementation of ISpectralMatch for writer tests.
        /// </summary>
        private class TestSpectralMatch : ISpectralMatch
        {
            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public int OneBasedScanNumber { get; }
            public double Score { get; }
            public double[]? QuantValues { get; set; }

            public TestSpectralMatch(string filePath, string fullSequence, string baseSequence, double score, int scanNumber)
            {
                FullFilePath = filePath;
                FullSequence = fullSequence;
                BaseSequence = baseSequence;
                Score = score;
                OneBasedScanNumber = scanNumber;
            }

            public int CompareTo(ISpectralMatch? other) => other == null ? 1 : Score.CompareTo(other.Score);
            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => Enumerable.Empty<IBioPolymerWithSetMods>();
        }

        #endregion
    }
}