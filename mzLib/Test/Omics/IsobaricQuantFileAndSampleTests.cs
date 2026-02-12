using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class IsobaricQuantFileAndSampleTests
    {
        /// <summary>
        /// Verify basic construction and that provided plex annotations are preserved.
        /// Also documents the canonical visible properties that callers read back from the instance.
        /// </summary>
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

        /// <summary>
        /// Ensure two instances with identical identity components (path, plex, fraction, techrep)
        /// compare equal and produce the same hash code.
        /// </summary>
        [Test]
        public void TwoObjects_WithSameIdentity_AreEqual()
        {
            var a = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>
            {
                new IsobaricQuantPlexAnnotation { Tag = "126", SampleName = "A" }
            });

            var b = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>
            {
                new IsobaricQuantPlexAnnotation { Tag = "127", SampleName = "B" }
            });

            Assert.That(a.Equals(b), Is.True, "Instances with identical identity fields should be equal.");
            Assert.That(b.Equals(a), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()), "Equal instances must have same hash code.");

            // operator overloads
            Assert.That(a == b, Is.True);
            Assert.That(a != b, Is.False);
        }

        /// <summary>
        /// When the identity differs by plex, fraction or technical replicate the instances must not be equal.
        /// </summary>
        [Test]
        public void TwoObjects_SamePathButDifferentMetadata_NotEqual()
        {
            var a = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var bDifferentPlex = new IsobaricQuantFileInfo("path/to/file.raw", "TMT11", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var bDifferentFraction = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 2, 1, new List<IsobaricQuantPlexAnnotation>());
            var bDifferentTechRep = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 1, 2, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(a.Equals(bDifferentPlex), Is.False);
            Assert.That(a.Equals(bDifferentFraction), Is.False);
            Assert.That(a.Equals(bDifferentTechRep), Is.False);

            Assert.That(a != bDifferentPlex, Is.True);
        }

        /// <summary>
        /// Instances with different file paths should not be equal (path is part of identity).
        /// </summary>
        [Test]
        public void Objects_WithDifferentFilePaths_AreNotEqual()
        {
            var a = new IsobaricQuantFileInfo("path/to/fileA.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var b = new IsobaricQuantFileInfo("path/to/fileB.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(a.Equals(b), Is.False, "Instances with different FullFilePathWithExtension should not be equal.");
            Assert.That(b.Equals(a), Is.False);
            Assert.That(a != b, Is.True);
        }

        /// <summary>
        /// Two equal objects must produce identical hash codes.
        /// </summary>
        [Test]
        public void GetHashCode_ReturnsSameValue_ForEqualObjects()
        {
            var a = new IsobaricQuantFileInfo("same/path/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());
            var b = new IsobaricQuantFileInfo("same/path/file.raw", "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(a.Equals(b), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
        }

        /// <summary>
        /// ToString should return the filename portion (no directory). Tests Windows-style path here.
        /// </summary>
        [Test]
        public void ToString_ReturnsFileNamePortion()
        {
            var filePath = Path.Combine("C:", "data", "experiments", "SampleQuantFile.raw");
            var fi = new IsobaricQuantFileInfo(filePath, "TMT10", 1, 1, new List<IsobaricQuantPlexAnnotation>());

            Assert.That(fi.ToString(), Is.EqualTo("SampleQuantFile.raw"));
        }

        // Edge cases

        /// <summary>
        /// Null and empty paths are normalized to empty string and treated equal.
        /// Also asserts that Annotations is never null even when passed null.
        /// </summary>
        [Test]
        public void NullAndEmptyPath_AreTreatedEquivalently()
        {
            var a = new IsobaricQuantFileInfo(null, "TMT10", 1, 1, null);
            var b = new IsobaricQuantFileInfo(string.Empty, "TMT10", 1, 1, null);

            Assert.That(a.FullFilePathWithExtension, Is.EqualTo(string.Empty));
            Assert.That(b.FullFilePathWithExtension, Is.EqualTo(string.Empty));
            Assert.That(a.Equals(b), Is.True);
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));

            Assert.That(a.Annotations, Is.Not.Null);
            Assert.That(b.Annotations, Is.Not.Null);
            Assert.That(a.Annotations.Count, Is.EqualTo(0));
            Assert.That(b.Annotations.Count, Is.EqualTo(0));
        }

        /// <summary>
        /// Constructor should normalize null plex and annotations inputs to non-null defaults.
        /// </summary>
        [Test]
        public void Constructor_NormalizesNullPlexAndAnnotations()
        {
            var fi = new IsobaricQuantFileInfo("p", null, 1, 1, null);

            Assert.That(fi.Plex, Is.EqualTo(string.Empty));
            Assert.That(fi.Annotations, Is.Not.Null);
            Assert.That(fi.Annotations.Count, Is.EqualTo(0));
        }

        /// <summary>
        /// Equals(object) with null should return false; typed IEquatable also should return false.
        /// </summary>
        [Test]
        public void Equals_WithNullObject_ReturnsFalse()
        {
            var a = new IsobaricQuantFileInfo("p", "TMT", 1, 1, null);

            Assert.That(a.Equals(null), Is.False);
            Assert.That(a == null, Is.False);
            Assert.That(null == a, Is.False);
        }

        /// <summary>
        /// Paths that differ only by separator style should be treated as equal because the implementation
        /// canonicalizes paths via Path.GetFullPath where possible.
        /// Also verify ToString still extracts the filename.
        /// </summary>
        [Test]
        public void PathSeparatorDifference_BecomesEqualAfterCanonicalization()
        {
            var forward = "dir/sub/file.raw";
            var back = @"dir\sub\file.raw";

            var f = new IsobaricQuantFileInfo(forward, "TMT", 1, 1, null);
            var b = new IsobaricQuantFileInfo(back, "TMT", 1, 1, null);

            // canonicalization with Path.GetFullPath should make these equivalent on most platforms
            Assert.That(f.Equals(b), Is.True, "Canonicalization should normalize separators and make the paths equal when they resolve to the same file.");
            Assert.That(f.ToString(), Is.EqualTo("file.raw"));
            Assert.That(b.ToString(), Is.EqualTo("file.raw"));
        }

        /// <summary>
        /// Paths that end with a directory separator yield an empty ToString() filename portion.
        /// Tests both forward and back slash behavior.
        /// </summary>
        [Test]
        public void PathEndingWithSlash_ToStringReturnsEmpty()
        {
            // Forward-slash ending
            var trailingForward = new IsobaricQuantFileInfo("dir/sub/", "TMT", 1, 1, null);
            Assert.That(trailingForward.ToString(), Is.EqualTo(string.Empty));

            // Backslash ending
            var trailingBack = new IsobaricQuantFileInfo(@"dir\sub\", "TMT", 1, 1, null);
            Assert.That(trailingBack.ToString(), Is.EqualTo(string.Empty));
        }

        /// <summary>
        /// If only a filename is provided, ToString returns that filename and equality behaves as expected.
        /// </summary>
        [Test]
        public void FilenameOnly_PathBehaves()
        {
            var fileOnly = new IsobaricQuantFileInfo("onlyfile.raw", "TMT", 1, 1, null);
            Assert.That(fileOnly.ToString(), Is.EqualTo("onlyfile.raw"));

            var fileOnly2 = new IsobaricQuantFileInfo("onlyfile.raw", "TMT", 1, 1, null);
            Assert.That(fileOnly.Equals(fileOnly2), Is.True);
            Assert.That(fileOnly == fileOnly2, Is.True);
        }

        /// <summary>
        /// Canonicalization does not change case; equality remains ordinal (case-sensitive).
        /// </summary>
        [Test]
        public void PathComparison_IsCaseSensitive()
        {
            var a = new IsobaricQuantFileInfo("PATH/To/File.RAW", "TMT10", 1, 1, null);
            var b = new IsobaricQuantFileInfo("path/to/file.raw", "TMT10", 1, 1, null);

            Assert.That(a.Equals(b), Is.False, "Ordinal comparison should treat different casing as not equal.");
        }

        /// <summary>
        /// Differences in fraction or technical replicate must make objects unequal.
        /// </summary>
        [Test]
        public void DifferentFractionOrTechnicalReplicate_AreNotEqual()
        {
            var baseInfo = new IsobaricQuantFileInfo("p.raw", "TMT10", 1, 1, null);
            var diffFraction = new IsobaricQuantFileInfo("p.raw", "TMT10", 2, 1, null);
            var diffTech = new IsobaricQuantFileInfo("p.raw", "TMT10", 1, 2, null);

            Assert.That(baseInfo.Equals(diffFraction), Is.False);
            Assert.That(baseInfo.Equals(diffTech), Is.False);
        }
    }
}
