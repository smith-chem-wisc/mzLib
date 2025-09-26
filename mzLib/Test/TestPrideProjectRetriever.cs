using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestPrideProjectRetriever
    {
        [Test]
        public void TestProjectRetrieve()
        {
            string projectAccession = "PXD048176";
            string outputFullFilePath = @"E:\junk\junk.txt";

            //UP000008595 is Uukuniemi virus (strain S23) (Uuk) which only has 4 proteins
            string returnedFilePath = PrideRetriever.RetrieveMassSpecProject("", outputFullFilePath);

            Assert.AreEqual(outputFullFilePath, returnedFilePath);
        }
        [Test]
        public void TestRetrieveProjectFileByName()
        {
            string projectAccession = "PXD048176";
            string filename = "d_atg1_d_atg11_proteome_data_analysis.7z";
            string outputDirectory = @"E:\junk";

            //UP000008595 is Uukuniemi virus (strain S23) (Uuk) which only has 4 proteins
            string returnedFilePath = PrideRetriever.RetrieveProjectFileByFilename(projectAccession,filename, outputDirectory);

            Assert.AreEqual(outputDirectory, returnedFilePath);
        }
        [Test]
        public void TestRetrieveProjectFileListByProjectAccession()
        {
            string projectAccession = "PXD048176";
            string outputDirectory = @"E:\junk";

            //UP000008595 is Uukuniemi virus (strain S23) (Uuk) which only has 4 proteins
            string returnedFilePath = PrideRetriever.RetrieveProjectFileListByProjectAccession(projectAccession,  outputDirectory);

            Assert.AreEqual(outputDirectory, returnedFilePath);
        }
        [Test]
        public void TestPrideFtp()
        {
            PrideRetriever.PrideFtp();
        }

        [Test]
        public void TestSimpleWebClientDownload()
        {
            string uri = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD048176/d_atg1_d_atg11_proteome_data_analysis.7z";
            string fullFilePath = @"E:\junk\PXD048176\d_atg1_d_atg11_proteome_data_analysis.7z";
            PrideRetriever.SimpleWebClientDownload(uri,fullFilePath);
        }

        [Test]
        public void DownloadUniProtProteomes()
        {
            var j = ProteinDbRetriever.DownloadAvailableUniProtProteomes(@"E:\junk");
        }
    }
}
