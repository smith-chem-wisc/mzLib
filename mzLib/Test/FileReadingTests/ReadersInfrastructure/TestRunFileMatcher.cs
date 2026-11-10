using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// <see cref="RunFileMatcher"/> called directly, as a reader without DIA-NN's records calls it. Each pairing rule
    /// is covered in depth through <see cref="DiaNnReportFile.FileNameToFilePath"/> in TestDiaNnResultFiles; these
    /// tests pin the shared entry point itself.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestRunFileMatcher
    {
        [Test]
        public void MatchAppliesEveryRuleOverPlainRunNames()
        {
            var runs = new[] { "S1", "sample_centroid", "A1", "R1_a", "R1_b" };
            var files = new List<string>
            {
                @"C:\data\s1.raw",          // exact, ignoring case and directory
                "/mnt/data/sample.mzML",     // the run carries a delimited suffix
                @"C:\data\A10.raw",          // "A1" must not match "A10"
                @"C:\data\R1.raw",           // both R1_a and R1_b extend it: contested, so neither keeps it
            };

            var matched = RunFileMatcher.Match(runs, files);

            Assert.That(matched, Has.Count.EqualTo(2));
            Assert.That(matched["S1"], Is.EqualTo(@"C:\data\s1.raw"));
            Assert.That(matched["sample_centroid"], Is.EqualTo("/mnt/data/sample.mzML"));
            Assert.That(matched.ContainsKey("A1"), Is.False);
            Assert.That(matched.ContainsKey("R1_a"), Is.False);
            Assert.That(matched.ContainsKey("R1_b"), Is.False);
        }

        [Test]
        public void MatchTakesEachRunNameOnce()
        {
            // A report repeats its run name on every row, and readers pass it in row by row
            var matched = RunFileMatcher.Match(new[] { "S1", "S1", "S1" }, new List<string> { @"C:\data\S1.raw" });

            Assert.That(matched, Has.Count.EqualTo(1));
            Assert.That(matched["S1"], Is.EqualTo(@"C:\data\S1.raw"));
        }

        [Test]
        [TestCase(@"C:\data\run1.raw", "run1")]
        [TestCase("/home/user/run1.mzML", "run1")]
        [TestCase("run1", "run1")]
        [TestCase("run1.wiff.scan", "run1.wiff")]
        [TestCase("", "")]
        [TestCase(null, "")]
        public void StripPathAndExtensionSplitsBothSeparators(string? path, string expected)
        {
            Assert.That(RunFileMatcher.StripPathAndExtension(path), Is.EqualTo(expected));
        }
    }
}
