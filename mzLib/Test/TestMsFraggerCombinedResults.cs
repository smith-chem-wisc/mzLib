using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;
using TopDownProteomics;
using OxyPlot;
using System.Diagnostics.CodeAnalysis;

namespace Test
{
    [ExcludeFromCodeCoverage]
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
        public void TestFileNameToFilePathWithParameter(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            List<string> fullFilePath = new List<string>();
            // these local files are not actually accessed, they are fillers to test the method
            string fullFilePath1 = @"E:\MadeleineH\Raw_Files\Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_11860_1x02nguL_8.raw";
            string fullFilePath2 = @"E:\MadeleineH\Raw_Files\Ex_AuLC1_30m_2D19_3_20um30cm_SPE50_15118120_OTOT_2215_HeYe_1.raw";
            fullFilePath.Add(fullFilePath1);
            fullFilePath.Add(fullFilePath2);

            List<string> results = ms.Results.Select(psm => psm.FileName).ToList();
            Dictionary<string, string> allFiles = ms.FileNameToFilePath(fullFilePath);
            List<string> filePaths = ms.ExperimentAnnotations.Select(psm => psm.File).ToList();

            foreach (var fileName in results)
            {
                Assert.That(allFiles.TryGetValue(fileName, out var output));
                Assert.That(filePaths.Contains(output));
            }
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

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults\experiment_annotation.tsv")]
        public void TestExperimentAnnotationFile(string path)
        {
            var testFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);

            ExperimentAnnotationFile experimentAnnotation = FileReader.ReadFile<ExperimentAnnotationFile>(testFilePath);

            experimentAnnotation.WriteResults(testFilePath);
            Assert.That(File.Exists(testFilePath));
        }
    }
}