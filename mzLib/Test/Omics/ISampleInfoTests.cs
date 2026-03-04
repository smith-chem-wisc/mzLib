using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.Omics
{
    /// <summary>
    /// Tests for ISampleInfo interface contract using real implementations.
    /// Note: Most ISampleInfo tests are in CrossClassComparisonTests.cs.
    /// This file only tests interface-level behavior not covered elsewhere.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class ISampleInfoTests
    {
        /// <summary>
        /// Verifies ISampleInfo implementations can be sorted in a collection.
        /// Critical: Sample ordering affects output file column order.
        /// </summary>
        [Test]
        public void ISampleInfo_Sorting_ProducesConsistentOrder()
        {
            var samples = new List<ISampleInfo>
            {
                new SpectraFileInfo(@"C:\b.raw", "Control", 2, 1, 0),
                new SpectraFileInfo(@"C:\a.raw", "Control", 1, 1, 0),
                new SpectraFileInfo(@"C:\a.raw", "Control", 1, 2, 0),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            // Should be ordered by: FilePath, then Condition, then BioRep, then TechRep, then Fraction
            Assert.That(samples[0].FullFilePathWithExtension, Is.EqualTo(@"C:\a.raw"));
            Assert.That(samples[0].TechnicalReplicate, Is.EqualTo(1));
            Assert.That(samples[1].TechnicalReplicate, Is.EqualTo(2));
            Assert.That(samples[2].FullFilePathWithExtension, Is.EqualTo(@"C:\b.raw"));
        }

        /// <summary>
        /// Verifies CompareTo handles null correctly.
        /// Critical: Prevents NullReferenceException during sorting.
        /// </summary>
        [Test]
        public void ISampleInfo_CompareTo_NullHandling()
        {
            var sample = new SpectraFileInfo(@"C:\test.raw", "Control", 1, 1, 0);

            // Non-null should come before null (negative return value)
            Assert.That(sample.CompareTo(null), Is.LessThan(0));
        }
    }
}