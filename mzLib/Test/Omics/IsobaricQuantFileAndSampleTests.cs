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
    }
}
