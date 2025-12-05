using MassSpectrometry;
using NUnit.Framework;
using Omics;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class IsobaricQuantFileAndSampleTests
    {
        [Test]
        public void TestIsobaricQuantFileAndPlexAnnotationCreation()
        {
            // Create annotations for the plex
            var annotations = new List<IsobaricQuantPlexAnnotation>
            {
                new IsobaricQuantPlexAnnotation
                {
                    Tag = "126",
                    SampleName = "Sample1",
                    Condition = "Control",
                    BiologicalReplicate = 1
                }
            };

            // Create an instance of IsobaricQuantFileInfo
            var fileInfo = new IsobaricQuantFileInfo(
                fullFilePathWithExtension: "SampleQuantFile.raw",
                plex: "TMT10",
                fraction: 1,
                technicalReplicate: 1,
                annotations: annotations);

            // Assert that the properties are set correctly
            Assert.That(fileInfo.FullFilePathWithExtension, Is.EqualTo("SampleQuantFile.raw"));
            Assert.That(fileInfo.Plex, Is.EqualTo("TMT10"));
            Assert.That(fileInfo.Fraction, Is.EqualTo(1));
            Assert.That(fileInfo.TechnicalReplicate, Is.EqualTo(1));
            Assert.That(fileInfo.Annotations, Is.Not.Null);
            Assert.That(fileInfo.Annotations.Count, Is.EqualTo(1));
            Assert.That(fileInfo.Annotations[0].Tag, Is.EqualTo("126"));
            Assert.That(fileInfo.Annotations[0].SampleName, Is.EqualTo("Sample1"));
            Assert.That(fileInfo.Annotations[0].Condition, Is.EqualTo("Control"));
            Assert.That(fileInfo.Annotations[0].BiologicalReplicate, Is.EqualTo(1));
        }

        [Test]
        public void TwoObjects_WithSameFilePath_AreEqual()
        {
            var a = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var b = new IsobaricQuantFileInfo("path/to/file.raw", "TMT11", 2, 2, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(a.Equals(b), Is.True, "Instances with identical FullFilePathWithExtension should be equal.");
            Assert.That(b.Equals(a), Is.True);
        }

        [Test]
        public void Objects_WithDifferentFilePaths_AreNotEqual()
        {
            var a = new IsobaricQuantFileInfo("path/to/fileA.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var b = new IsobaricQuantFileInfo("path/to/fileB.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(a.Equals(b), Is.False, "Instances with different FullFilePathWithExtension should not be equal.");
            Assert.That(b.Equals(a), Is.False);
        }

        [Test]
        public void GetHashCode_ReturnsSameValue_ForEqualObjects()
        {
            var a = new IsobaricQuantFileInfo("same/path/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var b = new IsobaricQuantFileInfo("same/path/file.raw", "TMTpro16", 2, 3, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(a.Equals(b), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        [Test]
        public void ToString_ReturnsFileNamePortion()
        {
            var filePath = @"C:\data\experiments\SampleQuantFile.raw";
            var fi = new IsobaricQuantFileInfo(filePath, "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(fi.ToString(), Is.EqualTo("SampleQuantFile.raw"));
        }
    }
}
