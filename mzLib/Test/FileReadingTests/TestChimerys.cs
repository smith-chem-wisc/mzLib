using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    public class TestChimerys
    {
        private static string directoryPath;

        [OneTimeSetUp]
        public void SetUp()
        {
            directoryPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ChimerysReadingWritingTests");
            Directory.CreateDirectory(directoryPath);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(directoryPath, true);
        }


        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_psms.tsv", 6)]
        public void TestChimerysPsmLoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysPsmFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_peptides.tsv", 5)]
        public void TestChimerysPeptidesLoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysPeptideFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }


        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_modified_peptides.tsv", 6)]
        public void TestChimerysModifiedPeptidesLoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysModifiedPeptideFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_precursors.tsv", 10)]
        public void TestChimerysPrecursorsLoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysPrecursorFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\Chimerysv4.3.0_MsaidPlatformv1.5.6_protein_groups.tsv", 4)]
        public void TestChimerysProteinGroupsLoadsAndCountCorrect(string path, int count)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            ChimerysProteinGroupFile file = new(filePath);
            Assert.That(file.Count(), Is.EqualTo(count));
            Assert.That(file.CanRead(path));
        }
    }
}
