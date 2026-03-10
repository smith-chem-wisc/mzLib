using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.Omics
{
    /// <summary>
    /// Tests for SpectraFileInfo class.
    /// Core equality and collection behavior is tested in CrossClassComparisonTests.cs.
    /// This file tests SpectraFileInfo-specific functionality.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SpectraFileInfoTests
    {
        /// <summary>
        /// Verifies all properties are correctly set by constructor.
        /// Critical: Ensures sample metadata is properly stored for quantification.
        /// </summary>
        [Test]
        public void Constructor_SetsAllProperties()
        {
            var sample = new SpectraFileInfo(
                @"C:\Data\Experiment\test_sample.raw", "Control", 1, 2, 3);

            Assert.Multiple(() =>
            {
                Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\Experiment\test_sample.raw"));
                Assert.That(sample.FilenameWithoutExtension, Is.EqualTo("test_sample"));
                Assert.That(sample.Condition, Is.EqualTo("Control"));
                Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
                Assert.That(sample.TechnicalReplicate, Is.EqualTo(2));
                Assert.That(sample.Fraction, Is.EqualTo(3));
            });
        }

        /// <summary>
        /// Verifies FilenameWithoutExtension correctly extracts filename from various path formats.
        /// Critical: Used for output file column headers and display.
        /// </summary>
        [Test]
        [TestCase(@"C:\Folder\MyFile.raw", "MyFile")]
        [TestCase(@"C:\Data\sample.test.raw", "sample.test")]
        [TestCase(@"C:\Data\sample.mzML", "sample")]
        [TestCase(@"C:\Data\samplefile", "samplefile")]
        [TestCase(@"\\server\share\data\sample.raw", "sample")]
        [TestCase("/home/user/data/sample.raw", "sample")]
        public void FilenameWithoutExtension_ExtractsCorrectly(string path, string expected)
        {
            var sample = new SpectraFileInfo(path, "Control", 1, 1, 0);
            Assert.That(sample.FilenameWithoutExtension, Is.EqualTo(expected));
        }

        /// <summary>
        /// Verifies ToString returns filename with extension (not full path).
        /// Critical: Used in output files and UI display.
        /// </summary>
        [Test]
        [TestCase(@"C:\Data\sample1.raw", "sample1.raw")]
        [TestCase(@"C:\Data\sample.mzML", "sample.mzML")]
        [TestCase(@"C:\Very\Deep\Path\file.raw", "file.raw")]
        public void ToString_ReturnsFilenameWithExtension(string path, string expected)
        {
            var sample = new SpectraFileInfo(path, "Control", 1, 1, 0);
            Assert.That(sample.ToString(), Is.EqualTo(expected));
        }

        /// <summary>
        /// Verifies CompareTo orders by Condition → BioRep → Fraction → TechRep → FilePath.
        /// Critical: Determines column order in quantification output files.
        /// </summary>
        [Test]
        public void CompareTo_OrdersByConditionThenBioRepThenFractionThenTechRepThenFilePath()
        {
            var samples = new List<SpectraFileInfo>
            {
                new(@"C:\z.raw", "Control", 2, 1, 0),
                new(@"C:\a.raw", "Control", 1, 2, 1),
                new(@"C:\b.raw", "Control", 1, 1, 0),
                new(@"C:\c.raw", "Alpha", 1, 1, 0),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            Assert.Multiple(() =>
            {
                // Alpha before Control (Condition first)
                Assert.That(samples[0].Condition, Is.EqualTo("Alpha"));

                // Within Control: BioRep 1 before BioRep 2
                Assert.That(samples[1].BiologicalReplicate, Is.EqualTo(1));
                Assert.That(samples[3].BiologicalReplicate, Is.EqualTo(2));

                // Within same BioRep: Fraction 0 before Fraction 1
                Assert.That(samples[1].Fraction, Is.EqualTo(0));
                Assert.That(samples[2].Fraction, Is.EqualTo(1));
            });
        }

        /// <summary>
        /// Verifies CompareTo handles null correctly.
        /// Critical: Prevents NullReferenceException during sorting.
        /// </summary>
        [Test]
        public void CompareTo_Null_ReturnsNegative()
        {
            var sample = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            Assert.That(sample.CompareTo(null), Is.LessThan(0));
        }

        /// <summary>
        /// Verifies identical samples compare as equal.
        /// Critical: Ensures sorting stability for duplicate entries.
        /// </summary>
        [Test]
        public void CompareTo_IdenticalValues_ReturnsZero()
        {
            var sample1 = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);

            Assert.That(sample1.CompareTo(sample2), Is.EqualTo(0));
        }

        /// <summary>
        /// Verifies SpectraFileInfo deduplicates correctly in HashSet.
        /// Critical: Prevents duplicate sample entries in quantification.
        /// </summary>
        [Test]
        public void HashSet_DeduplicatesCorrectly()
        {
            var sample1 = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var duplicate = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var different = new SpectraFileInfo(@"C:\other.raw", "Control", 1, 1, 0);

            var set = new HashSet<SpectraFileInfo> { sample1, duplicate, different };

            Assert.That(set.Count, Is.EqualTo(2));
        }

        /// <summary>
        /// Verifies SpectraFileInfo works correctly as Dictionary key.
        /// Critical: IntensitiesBySample uses ISampleInfo as key.
        /// </summary>
        [Test]
        public void Dictionary_WorksAsKey()
        {
            var sample1 = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);
            var duplicate = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);

            var dict = new Dictionary<SpectraFileInfo, double> { { sample1, 1000.0 } };

            Assert.That(dict.TryGetValue(duplicate, out var value), Is.True);
            Assert.That(value, Is.EqualTo(1000.0));
        }
    }
}