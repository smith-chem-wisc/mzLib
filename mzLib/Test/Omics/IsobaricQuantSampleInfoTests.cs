using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics
{
    /// <summary>
    /// Tests for IsobaricQuantSampleInfo class.
    /// Tests unique behavior not covered in CrossClassComparisonTests.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class IsobaricQuantSampleInfoTests
    {
        /// <summary>
        /// Verifies all properties are correctly set by constructor.
        /// Critical: Ensures TMT/iTRAQ channel metadata is properly stored.
        /// </summary>
        [Test]
        public void Constructor_SetsAllProperties()
        {
            var sample = new IsobaricQuantSampleInfo(
                @"C:\Data\test.raw", "Control", 1, 2, 3, 4, "127N", 127.124761, true);

            Assert.Multiple(() =>
            {
                Assert.That(sample.FullFilePathWithExtension, Is.EqualTo(@"C:\Data\test.raw"));
                Assert.That(sample.Condition, Is.EqualTo("Control"));
                Assert.That(sample.BiologicalReplicate, Is.EqualTo(1));
                Assert.That(sample.TechnicalReplicate, Is.EqualTo(2));
                Assert.That(sample.Fraction, Is.EqualTo(3));
                Assert.That(sample.PlexId, Is.EqualTo(4));
                Assert.That(sample.ChannelLabel, Is.EqualTo("127N"));
                Assert.That(sample.ReporterIonMz, Is.EqualTo(127.124761));
                Assert.That(sample.IsReferenceChannel, Is.True);
            });
        }

        /// <summary>
        /// Verifies UniqueIdentifier is computed from FilePath and ChannelLabel only.
        /// Critical: UniqueIdentifier is used as HashCode; determines collection behavior.
        /// </summary>
        [Test]
        public void UniqueIdentifier_BasedOnFilePathAndChannelLabel()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sameIdentity = new IsobaricQuantSampleInfo(@"C:\a.raw", "B", 99, 99, 99, 99, "126", 999.0, true);
            var differentChannel = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "127N", 126.0, false);
            var differentFile = new IsobaricQuantSampleInfo(@"C:\b.raw", "A", 1, 1, 0, 1, "126", 126.0, false);

            Assert.Multiple(() =>
            {
                Assert.That(sample1.UniqueIdentifier, Is.EqualTo(sameIdentity.UniqueIdentifier));
                Assert.That(sample1.UniqueIdentifier, Is.Not.EqualTo(differentChannel.UniqueIdentifier));
                Assert.That(sample1.UniqueIdentifier, Is.Not.EqualTo(differentFile.UniqueIdentifier));
                Assert.That(sample1.GetHashCode(), Is.EqualTo(sample1.UniqueIdentifier));
            });
        }

        /// <summary>
        /// Verifies ToString returns filename (without extension) and channel label.
        /// Critical: Used in output file headers for intensity columns.
        /// </summary>
        [Test]
        public void ToString_ReturnsFileNameAndChannelLabel()
        {
            var sample = new IsobaricQuantSampleInfo(@"C:\Path\To\Experiment.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            Assert.That(sample.ToString(), Is.EqualTo("Experiment_126"));
        }

        /// <summary>
        /// Verifies CompareTo orders by Condition → BioRep → Fraction → TechRep → FilePath.
        /// Critical: Determines column order in quantification output files.
        /// </summary>
        [Test]
        public void CompareTo_OrdersByConditionThenBioRepThenFractionThenTechRepThenFilePath()
        {
            var samples = new List<IsobaricQuantSampleInfo>
            {
                new(@"C:\z.raw", "Beta", 1, 1, 0, 1, "126", 126.0, false),
                new(@"C:\a.raw", "Alpha", 2, 1, 0, 1, "127N", 127.0, false),
                new(@"C:\b.raw", "Alpha", 1, 2, 1, 1, "128C", 128.0, false),
                new(@"C:\c.raw", "Alpha", 1, 1, 0, 1, "129N", 129.0, false),
            };

            samples.Sort((a, b) => a.CompareTo(b));

            Assert.Multiple(() =>
            {
                // Alpha before Beta (Condition)
                Assert.That(samples[0].Condition, Is.EqualTo("Alpha"));
                Assert.That(samples[3].Condition, Is.EqualTo("Beta"));

                // Within Alpha: BioRep 1 before BioRep 2
                Assert.That(samples[0].BiologicalReplicate, Is.EqualTo(1));
                Assert.That(samples[2].BiologicalReplicate, Is.EqualTo(2));

                // Within same BioRep: Fraction 0 before Fraction 1
                Assert.That(samples[0].Fraction, Is.EqualTo(0));
                Assert.That(samples[1].Fraction, Is.EqualTo(1));
            });
        }

        /// <summary>
        /// Verifies CompareTo returns 0 for samples with same FilePath and ChannelLabel.
        /// Critical: Consistent with Equals behavior for sorting stability.
        /// </summary>
        [Test]
        public void CompareTo_SameIdentity_ReturnsZero()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "Treatment", 99, 99, 99, 99, "126", 999.0, true);

            Assert.That(sample1.CompareTo(sample2), Is.EqualTo(0));
        }

        /// <summary>
        /// Verifies all standard TMT channel labels work correctly.
        /// Critical: Ensures real-world TMT experiments are properly handled.
        /// </summary>
        [Test]
        public void TMTChannels_AllDistinct()
        {
            var channels = new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" };
            var samples = channels.Select((ch, i) =>
                new IsobaricQuantSampleInfo(@"C:\test.raw", "Control", 1, 1, 0, 1, ch, 126.0 + i * 0.5, false)).ToList();

            var set = new HashSet<IsobaricQuantSampleInfo>(samples);

            Assert.That(set.Count, Is.EqualTo(11), "All 11 TMT channels should be distinct");
        }
    }
}