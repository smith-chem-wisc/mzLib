using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;
using NUnit.Framework;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using FlashLFQ;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestPrideProjectRetriever
    {
        [Test]
        public void TestProjectRetrieve()
        {
            var outputFullFileFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            var outputDir = Path.Combine(outputFullFileFolder, @"PrideRetrieveTest");
            Directory.CreateDirectory(outputDir);
            string projectAccession = "PXD047650";
            PrideRetriever newPride = new PrideRetriever(projectAccession, outputDir);
            List<PRIDEEntry> pRIDEEntry = newPride.GetPrideEntries();
            Assert.AreEqual(3, pRIDEEntry.Count);
            Directory.Delete(outputDir, true);
        }

        [Test]
        public void TestRetrieveProjectFileByName()
        {
            var outputFullFileFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            var outputDir = Path.Combine(outputFullFileFolder, @"PrideRetrieveTest");
            Directory.CreateDirectory(outputDir);
            string projectAccession = "PXD047650";

            PrideRetriever newPrideProject = new PrideRetriever(projectAccession, outputDir);
            string returnedFilePath = newPrideProject.RetrieveProjectFileByFilename(newPrideProject.GetPrideFileNames().Last(), outputDir);
            String outputDirectory = Path.Combine(outputDir, newPrideProject.GetPrideFileNames().Last());

            Assert.AreEqual(outputDirectory, returnedFilePath);
            Directory.Delete(outputDir, true);
        }

        [Test]
        public void TestRetrieveExperimentalInformation()
        {
            var outputFullFileFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            var outputDir = Path.Combine(outputFullFileFolder, @"PrideRetrieveTest");
            Directory.CreateDirectory(outputDir);
            string projectAccession = "PXD047650";

            PrideRetriever newPrideProject = new PrideRetriever(projectAccession, outputDir);
            PrideProject experiment = newPrideProject.GetPrideProjectInformation();

            Assert.AreEqual(projectAccession, experiment.Accession);
            Directory.Delete(outputDir, true);
        }

        [Test]
        public void TestPrideFtp()
        {
            PrideRetriever.PrideFtp();
        }

        [Test]
        public void TestSimpleWebClientDownload()
        {
            var outputFullFileFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            var outputDir = Path.Combine(outputFullFileFolder, @"PrideRetrieveTest");
            Directory.CreateDirectory(outputDir);
            string projectAccession = "PXD047650";

            PrideRetriever newPrideProject = new PrideRetriever(projectAccession, outputDir);
            String testFileName = newPrideProject.GetPrideFileNames().Last();
            String uri = newPrideProject.GetFtpLinkByFileNmae(testFileName);

            string fullFilePath = Path.Combine(outputDir, testFileName);
            string expectedFile = Path.Combine(outputDir, @"proteinGroups.xlsx");
            PrideRetriever.SimpleWebClientDownload(uri, fullFilePath);
            Assert.AreEqual(fullFilePath, expectedFile);
            Directory.Delete(outputDir, true);
        }

        [Test]
        public void TestGetAllFilesByAccession()
        {
            var outputFullFileFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData");
            var outputDir = Path.Combine(outputFullFileFolder, @"PrideRetrieveTest");
            Directory.CreateDirectory(outputDir);
            string projectAccession = "PXD047650";

            PrideRetriever newPrideProject = new PrideRetriever(projectAccession, outputDir);
            List<String> testFileName = newPrideProject.GetAllFilesByAccession(projectAccession, outputDir);

            Assert.AreEqual(testFileName.Count, 3);
            Directory.Delete(outputDir, true);
        }


        [Test]
        public void TestSearchByKeywordsAndFilters()
        {
            int pageSize = 10;
            var accessionList = PrideRetriever.SearchByKeywordsAndFilters("human", "", pageSize, 0, "", "ASC", "submission_date");
            Assert.AreEqual(accessionList.Count, pageSize);
        }

    }
}
