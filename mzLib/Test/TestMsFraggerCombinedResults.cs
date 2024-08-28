using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;
using TopDownProteomics;
using OxyPlot;

namespace Test
{
    internal class TestMsFraggerCombinedResults
    {
        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults")]
        public void TestLoadResultsCount(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            Assert.That(ms.AllPsmFiles.Count.Equals(2));
            Assert.That(ms.Results.Count.Equals(8));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults")]
        public void TestLoadResults(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            List<string> results = ms.Results.Select(psm => psm.FileName).ToList();

            Assert.That((results.Count(s => s.Contains("A_1"))).Equals(4));
            Assert.That((results.Count(s => s.Contains("A_2"))).Equals(4));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults")]
        public void TestFileNameToFilePathWithoutParameter(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            List<string> results = ms.Results.Select(psm => psm.FileName).ToList();
            Dictionary<string, string> allFiles = ms.FileNameToFilePath();
            List<string> filePaths = ms.ExperimentAnnotations.Select(psm => psm.File).ToList();

            foreach (var fileName in results)
            {
                Assert.That(allFiles.TryGetValue(fileName, out var output));
                Assert.That(filePaths.Contains(output));
            }
        }
    }
}