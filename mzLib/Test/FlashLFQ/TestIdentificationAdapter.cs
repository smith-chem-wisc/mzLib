using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;

namespace Test.FlashLFQ
{
    public class TestIdentificationAdapter
    {
        //[Test]
        //public void LocalFileTest()
        //{
        //    List<SpectraFileInfo> spectraFiles = new List<SpectraFileInfo>
        //    {
        //        new SpectraFileInfo(@"D:\FlashLFQVignette-selected\09-04-18_EcoliSpikeInSingleShot1x.raw",
        //            "1x", 0, 0, 0),
        //        new SpectraFileInfo(@"D:\FlashLFQVignette-selected\09-04-18_EcoliSpikeInSingleShot2x.raw",
        //            "2x", 0, 0, 0),
        //    };
        //    string pepPath = @"D:\FlashLFQVignette-selected\MMv1p1p0_Search\Task1-SearchTask\AllPeptides.psmtsv";
        //    string psmPath = @"D:\FlashLFQVignette-selected\MMv1p1p0_Search\Task1-SearchTask\AllPSMs.psmtsv";

        //    IQuantifiableResultFile psmFile = FileReader.ReadQuantifiableResultFile(pepPath);
        //    var ids = psmFile.MakeIdentifications(spectraFiles, usePepQValue: true);
        //    int placeholder = 0;
        //}
        

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv")]
        public void TestAddProteinGroupInfoCorrect(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            List<Identification> identifications = new List<Identification>();
            List<SpectraFileInfo> spectraFiles = new List<SpectraFileInfo>
                {
                    new SpectraFileInfo(@"D:\Projects\Chimeras\Mann_11cell_analysis\Hela\MsFragger\Hela_1_1\interact-20100611_Velos1_TaGe_SA_Hela_1.pep.xml", 
                                        "", 0, 0, 0),
                };
            identifications = MzLibExtensions.MakeIdentifications(file, spectraFiles);

            // list should contain five elements
            Assert.That(identifications.Count, Is.EqualTo(5));
            // one protein associated with given results, list should only contain this one element 
            Assert.That(identifications[0].ProteinGroups.Count, Is.EqualTo(1));
            // two proteins associated with given results, list should contain two elements
            Assert.That(identifications[2].ProteinGroups.Count, Is.EqualTo(2));

            Identification identification1 = identifications[0];
            Assert.That(identification1.BaseSequence, Is.EqualTo("KPVGAAK"));
            Assert.That(identification1.ModifiedSequence, Is.EqualTo("KPVGAAK"));
            Assert.That(identification1.Ms2RetentionTimeInMinutes, Is.EqualTo(1.9398));
            Assert.That(identification1.MonoisotopicMass, Is.EqualTo(669.4173));
            Assert.That(identification1.PrecursorChargeState, Is.EqualTo(2));

            HashSet<ProteinGroup> proteinGroups = identification1.ProteinGroups;
            ProteinGroup proteinGroup1 = proteinGroups.First();
            Assert.That(proteinGroup1.ProteinGroupName, Is.EqualTo("P16403"));
            Assert.That(proteinGroup1.GeneName, Is.EqualTo("H12"));
            Assert.That(proteinGroup1.Organism, Is.EqualTo("HUMAN"));

            Identification identification5 = identifications[4];
            Assert.That(identification5.BaseSequence, Is.EqualTo("VVTHGGR"));
            Assert.That(identification5.ModifiedSequence, Is.EqualTo("VVTHGGR"));
            Assert.That(identification5.Ms2RetentionTimeInMinutes, Is.EqualTo(19.114));
            Assert.That(identification5.MonoisotopicMass, Is.EqualTo(724.398));
            Assert.That(identification5.PrecursorChargeState, Is.EqualTo(2));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv")]
        public void TestFileNametoFilePath(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            string fileName = file.First().FileName;

            List<string> fullFilePath = new List<string>();
            string fullFilePath1 = @"D:\Projects\Chimeras\Mann_11cell_analysis\RawData\interact-20100611_Velos1_TaGe_SA_Hela_1.raw";
            string fullFilePath2 = @"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv";
            fullFilePath.Add(fullFilePath1);
            fullFilePath.Add(fullFilePath2);

            Dictionary<string, string> allFiles = file.FileNameToFilePath(fullFilePath);

            Assert.That(allFiles.TryGetValue(fileName, out var output));
            Assert.AreEqual(output, fullFilePath1);
            Assert.That(!allFiles.ContainsValue(fullFilePath2));
        }

        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\SmallCalibratibleYeastFragger_psm.tsv")]
        public void TestFileNametoFilePathLocalPath(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);
            string fileName = file.First().FileName;

            List<string> fullFilePath = new List<string>();
            string rawFilePath = @"DataFiles\SmallCalibratibleYeast.mzml";
            fullFilePath.Add(rawFilePath);

            Dictionary<string, string> allFiles = file.FileNameToFilePath(fullFilePath);

            Assert.That(allFiles.TryGetValue(fileName, out var output));
            Assert.AreEqual(output, rawFilePath);
        }

        [Test]
        public void MakeIdentifications_ShouldReturnEmptyList_WhenNoQuantifiableRecords()
        {
            // Arrange
            var quantifiable = new MockQuantifiableResultFile(new List<IQuantifiableRecord>());
            var spectraFiles = new List<SpectraFileInfo>();

            // Act
            var result = MzLibExtensions.MakeIdentifications(quantifiable, spectraFiles);

            // Assert
            Assert.IsEmpty(result);
        }

        [Test]
        public void MakeIdentifications_ShouldReturnIdentifications_WhenQuantifiableRecordsExist()
        {
            // Arrange
            var quantifiableRecords = new List<IQuantifiableRecord>
                {
                    new MockQuantifiableRecord
                    {
                        BaseSequence = "BASESEQ",
                        FullSequence = "MODSEQ",
                        RetentionTime = 5.0,
                        MonoisotopicMass = 500.0,
                        ChargeState = 2,
                        FileName = "file1.mzML",
                        ProteinGroupInfos = new List<(string proteinAccessions, string geneName, string organism)>
                        {
                            ("P1", "Gene1", "Organism1")
                        }
                    },
                    new MockQuantifiableRecord
                    {
                        BaseSequence = "BASESEQ2",
                        FullSequence = "MODSEQ2",
                        RetentionTime = 10.0,
                        MonoisotopicMass = 1000.0,
                        ChargeState = 3,
                        FileName = "file2.mzML",
                        ProteinGroupInfos = new List<(string proteinAccessions, string geneName, string organism)>
                        {
                            ("P2", "Gene2", "Organism2")
                        }
                    }
                };
            var quantifiable = new MockQuantifiableResultFile(quantifiableRecords);
            var spectraFiles = new List<SpectraFileInfo>
                {
                    new SpectraFileInfo("file1.mzML", "", 0, 0, 0),
                    new SpectraFileInfo("file2.mzML", "", 0, 0, 0)
                };

            // Act
            var result = MzLibExtensions.MakeIdentifications(quantifiable, spectraFiles);

            // Assert
            Assert.AreEqual(2, result.Count);
            var identification1 = result[0];
            Assert.AreEqual("BASESEQ", identification1.BaseSequence);
            Assert.AreEqual("MODSEQ", identification1.ModifiedSequence);
            Assert.AreEqual(5.0, identification1.Ms2RetentionTimeInMinutes);
            Assert.AreEqual(500.0, identification1.MonoisotopicMass);
            Assert.AreEqual(2, identification1.PrecursorChargeState);
            Assert.AreEqual(1, identification1.ProteinGroups.Count);
            Assert.AreEqual("P1", identification1.ProteinGroups.First().ProteinGroupName);
            Assert.AreEqual("Gene1", identification1.ProteinGroups.First().GeneName);
            Assert.AreEqual("Organism1", identification1.ProteinGroups.First().Organism);
            Assert.AreEqual(spectraFiles[0], identification1.FileInfo);

            var identification2 = result[1];
            Assert.AreEqual("BASESEQ2", identification2.BaseSequence);
            Assert.AreEqual("MODSEQ2", identification2.ModifiedSequence);
            Assert.AreEqual(10.0, identification2.Ms2RetentionTimeInMinutes);
            Assert.AreEqual(1000.0, identification2.MonoisotopicMass);
            Assert.AreEqual(3, identification2.PrecursorChargeState);
            Assert.AreEqual(1, identification2.ProteinGroups.Count);
            Assert.AreEqual("P2", identification2.ProteinGroups.First().ProteinGroupName);
            Assert.AreEqual("Gene2", identification2.ProteinGroups.First().GeneName);
            Assert.AreEqual("Organism2", identification2.ProteinGroups.First().Organism);
            Assert.AreEqual(spectraFiles[1], identification2.FileInfo);
        }

        [Test]
        public void SpectraFileNotFound()
        {
            var quantifiableRecords = new List<IQuantifiableRecord>
                {
                    new MockQuantifiableRecord
                    {
                        BaseSequence = "BASESEQ",
                        FullSequence = "MODSEQ",
                        RetentionTime = 5.0,
                        MonoisotopicMass = 500.0,
                        ChargeState = 2,
                        FileName = "file1.mzML",
                        ProteinGroupInfos = new List<(string proteinAccessions, string geneName, string organism)>
                        {
                            ("P1", "Gene1", "Organism1")
                        }
                    }
                };
            var quantifiable = new MockQuantifiableResultFile(quantifiableRecords);
            var spectraFiles = new List<SpectraFileInfo>
                {
                    new SpectraFileInfo("file2.mzML", "", 0, 0, 0)
                };

            try
            {
                var result = MzLibExtensions.MakeIdentifications(quantifiable, spectraFiles);
            }
            catch (Exception ex)
            {
                Assert.AreEqual("Spectra file not found for file name: file1.mzML", ex.Message);
            }
        }

        [Test]
        public void TestMakeIdentificationsWithLegacyPsmTsvInput()
        {
            string psmFilename = "AllPSMs2.psmtsv";

            var myDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");
            var pathOfIdentificationFile = Path.Combine(myDirectory, psmFilename);
            var pathOfMzml = Path.Combine(myDirectory, "SmallCalibratible_Yeast.mzML");
            SpectraFileInfo sfi = new SpectraFileInfo(pathOfMzml, "A", 1, 1, 1);

            PsmFromTsvFile results = new PsmFromTsvFile(pathOfIdentificationFile);
            Assert.That(results.Results.Count, Is.EqualTo(89));

            IQuantifiableResultFile quantifiableResultFile = FileReader.ReadQuantifiableResultFile(pathOfIdentificationFile);
            List<Identification> ids = MzLibExtensions.MakeIdentifications(quantifiableResultFile, new List<SpectraFileInfo> { sfi });
            Assert.That(ids.Count, Is.EqualTo(89));

            var testResult = (PsmFromTsv)quantifiableResultFile.GetQuantifiableResults().First();

            Assert.That(testResult != null);
            Assert.That(Double.IsNaN(testResult.PEP));
            Assert.That(Double.IsNaN(testResult.PEP_QValue));
            Assert.That(!Double.IsNaN(testResult.RetentionTime));
        }

        [Test]
        public static void MakeIdentificationsFromNewPsmtsv()
        {
            string psmFilename = "AllPSMsWNewFormat.psmtsv";

            var myDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");
            var pathOfIdentificationFile = Path.Combine(myDirectory, psmFilename);
            var pathOfMzml = Path.Combine(myDirectory, "SmallCalibratibleYeast.mzML");
            SpectraFileInfo sfi = new SpectraFileInfo(pathOfMzml, "A", 1, 1, 1);

            PsmFromTsvFile results = new PsmFromTsvFile(pathOfIdentificationFile);
            Assert.That(results.Results.Count, Is.EqualTo(109));

            IQuantifiableResultFile quantifiableResultFile = FileReader.ReadQuantifiableResultFile(pathOfIdentificationFile);
            List<Identification> ids = quantifiableResultFile.MakeIdentifications(new List<SpectraFileInfo> { sfi });
            Assert.That(ids.Count, Is.EqualTo(109));

            var decoyId = ids.FirstOrDefault(r => r.ModifiedSequence == "GAAVVQK");
            Assert.IsTrue(decoyId.IsDecoy);
            Assert.IsFalse(decoyId.UseForProteinQuant);
            Assert.That(decoyId.QValue, Is.EqualTo(0.013889).Within(0.000001));
            Assert.That(decoyId.PsmScore, Is.EqualTo(6.218).Within(0.001));


            ids = quantifiableResultFile.MakeIdentifications(new List<SpectraFileInfo> { sfi }, usePepQValue: true);
            Assert.That(ids.Count, Is.EqualTo(109));
            decoyId = ids.FirstOrDefault(r => r.ModifiedSequence == "GAAVVQK");
            Assert.IsTrue(decoyId.IsDecoy);
            Assert.IsFalse(decoyId.UseForProteinQuant);
            Assert.That(decoyId.QValue, Is.EqualTo(2).Within(0.000001));
            Assert.That(decoyId.PsmScore, Is.EqualTo(6.218).Within(0.001));
        }
    }

    // Mock classes for testing
    public class MockQuantifiableResultFile : IResultFile, IQuantifiableResultFile
    {
        private readonly List<IQuantifiableRecord> _quantifiableRecords;

        public string FilePath { get; set; }

        public MockQuantifiableResultFile(List<IQuantifiableRecord> quantifiableRecords)
        {
            _quantifiableRecords = quantifiableRecords;
        }

        public IEnumerable<IQuantifiableRecord> GetQuantifiableResults()
        {
            return _quantifiableRecords;
        }

        public Dictionary<string, string> FileNameToFilePath(List<string> fullFilePaths)
        {
            var dict = new Dictionary<string, string>();
            foreach (var path in fullFilePaths)
            {
                dict[path] = path;
            }
            return dict;
        }

        // public string FilePath => throw new System.NotImplementedException();
        public SupportedFileType FileType => SupportedFileType.Mgf;
        public Software Software { get => Software.Crux; set => throw new System.NotImplementedException(); }
        //public string FilePath { get => null; set => throw new System.NotImplementedException(); }

        public void LoadResults() => throw new System.NotImplementedException();
        public void WriteResults(string path) => throw new System.NotImplementedException();
    }

    public class MockQuantifiableRecord : IQuantifiableRecord
    {
        public string BaseSequence { get; set; }
        public string FullSequence { get; set; }
        public double RetentionTime { get; set; }
        public double MonoisotopicMass { get; set; }
        public int ChargeState { get; set; }
        public List<(string proteinAccessions, string geneName, string organism)> ProteinGroupInfos { get; set; }
        public string FileName { get; set; }
        public bool IsDecoy { get; set; }
    }

}