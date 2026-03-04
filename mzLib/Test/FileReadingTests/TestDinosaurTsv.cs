using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.FileReadingTests
{
    internal class TestDinosaurTsv
    {

        private static string directoryPath;
        private static string filePath => Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\DinoSnippet.features.tsv");

        [OneTimeSetUp]
        public void SetUp()
        {
            directoryPath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ReadingWritingTests");
            Directory.CreateDirectory(directoryPath);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(directoryPath, true);
        }

        [Test]
        public void TestDinosaurFeatureProperties()
        {
            var feature = new DinosaurFeature
            {
                Mz = 123.45,
                MostAbundantMz = 123.46,
                Charge = 2,
                RtStart = 10.1,
                RtApex = 10.5,
                RtEnd = 10.9,
                Fwhm = 0.2,
                NIsotopes = 3,
                NScans = 5,
                AveragineCorr = 0.98,
                Mass = 246.9,
                MassCalib = 247.0,
                IntensityApex = 10000.0,
                IntensitySum = 50000.0
            };

            Assert.That(feature.Mz, Is.EqualTo(123.45));
            Assert.That(feature.MostAbundantMz, Is.EqualTo(123.46));
            Assert.That(feature.Charge, Is.EqualTo(2));
            Assert.That(feature.RtStart, Is.EqualTo(10.1));
            Assert.That(feature.RtApex, Is.EqualTo(10.5));
            Assert.That(feature.RtEnd, Is.EqualTo(10.9));
            Assert.That(feature.Fwhm, Is.EqualTo(0.2));
            Assert.That(feature.NIsotopes, Is.EqualTo(3));
            Assert.That(feature.NScans, Is.EqualTo(5));
            Assert.That(feature.AveragineCorr, Is.EqualTo(0.98));
            Assert.That(feature.Mass, Is.EqualTo(246.9));
            Assert.That(feature.MassCalib, Is.EqualTo(247.0));
            Assert.That(feature.IntensityApex, Is.EqualTo(10000.0));
            Assert.That(feature.IntensitySum, Is.EqualTo(50000.0));
        }

        [Test]
        public static void TestParameterlessConstructor()
        {
            DinosaurTsvFile x = new DinosaurTsvFile();
            Assert.That(x.FileType, Is.EqualTo(SupportedFileType.Tsv_Dinosaur), "File type should be Tsv_Dinosaur.");
        }

        [Test]
        public void TestDinosaurTsvFileLoadResults()
        {
            var dinoFile = new DinosaurTsvFile(filePath);
            dinoFile.LoadResults();
            Assert.That(dinoFile.Results, Is.Not.Null);
            Assert.That(dinoFile.Results.Count, Is.EqualTo(9), "Expected 9 entries in the Dinosaur .feature.tsv file.");

            // Check that all required properties are populated for the first entry
            var first = dinoFile.Results.First();
            Assert.That(first.Mz, Is.EqualTo(440.2366771));
            Assert.That(first.MostAbundantMz, Is.EqualTo(440.2366771));
            Assert.That(first.Charge, Is.EqualTo(3));
            Assert.That(first.RtStart, Is.EqualTo(32.67795675));
            Assert.That(first.RtApex, Is.EqualTo(32.72369833));
            Assert.That(first.RtEnd, Is.EqualTo(32.73719507));
            Assert.That(first.Fwhm, Is.EqualTo(1.021054854));
            Assert.That(first.NIsotopes, Is.EqualTo(2));
            Assert.That(first.NScans, Is.EqualTo(5));
            Assert.That(first.AveragineCorr, Is.EqualTo(0.992839779));
            Assert.That(first.Mass, Is.EqualTo(1317.687838));
            Assert.That(first.MassCalib, Is.EqualTo(1317.687788));
            Assert.That(first.IntensityApex, Is.EqualTo(8809.755127));
            Assert.That(first.IntensitySum, Is.EqualTo(52349.81689));
        }

        [Test]
        public void TestDinosaurIMs1FeatureFileImplementation()
        {
            var dinoFile = new DinosaurTsvFile(filePath);
            dinoFile.LoadResults();

            // Check that all required properties are populated for the first entry
            var first = dinoFile.Results.First();

            var firstIMs1Feature = dinoFile.GetMs1Features().First();
            Assert.That(first, Is.EqualTo(firstIMs1Feature), "First entry from Results should match first from GetMs1Features().");
        }

        [Test]
        public void TestDinosaurTsvFileWriteResults()
        {
            var dinoFile = new DinosaurTsvFile(filePath);
            dinoFile.LoadResults();
            string outPath = Path.Combine(directoryPath, "DinoSnippet.feature.tsv");
            dinoFile.WriteResults(outPath);

            Assert.That(File.Exists(outPath), Is.True, "Output file was not created.");

            // Reload and check count matches
            var dinoFileReloaded = new DinosaurTsvFile(outPath);
            dinoFileReloaded.LoadResults();
            Assert.That(dinoFileReloaded.Results.Count, Is.EqualTo(9), "Reloaded file should have 9 entries.");

            outPath = Path.Combine(directoryPath, "DinoSnippetNoExtension");
            dinoFile.WriteResults(outPath);

            Assert.That(File.Exists(outPath + ".feature.tsv"), Is.True, "Output file was not created.");
        }
    }
}
