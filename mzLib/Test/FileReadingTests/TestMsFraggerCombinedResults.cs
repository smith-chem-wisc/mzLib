using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using TopDownProteomics;
using System.Diagnostics.CodeAnalysis;

namespace Test.FileReadingTests
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

            NUnit.Framework.Assert.That(ms.AllPsmFiles.Count.Equals(2));
            NUnit.Framework.Assert.That(ms.Results.Count.Equals(8));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults")]
        public void TestLoadResults(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerCombinedResults ms = new MsFraggerCombinedResults(filePath);
            ms.LoadResults();

            List<string> results = ms.Results.Select(psm => psm.FileName).ToList();

            NUnit.Framework.Assert.That(results.Count(s => s.Contains("A_1")).Equals(4));
            NUnit.Framework.Assert.That(results.Count(s => s.Contains("A_2")).Equals(4));
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
                NUnit.Framework.Assert.That(allFiles.TryGetValue(fileName, out var output));
                NUnit.Framework.Assert.That(filePaths.Contains(output));
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
                NUnit.Framework.Assert.That(allFiles.TryGetValue(fileName, out var output));
                NUnit.Framework.Assert.That(filePaths.Contains(output));
            }
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults\experiment_annotation.tsv")]
        public void TestExperimentAnnotationFile(string path)
        {
            string fileToWrite = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\EditedMSFraggerResults\copy_experiment_annotation.tsv");
            if (File.Exists(fileToWrite))
            {
                File.Delete(fileToWrite);
            }

            string fileToRead = Path.Combine(TestContext.CurrentContext.TestDirectory, path);

            ExperimentAnnotationFile experimentAnnotation = FileReader.ReadFile<ExperimentAnnotationFile>(fileToRead);

            experimentAnnotation.WriteResults(fileToWrite);
            NUnit.Framework.Assert.That(File.Exists(fileToWrite));

            File.Delete(fileToWrite);
        }
    }
}