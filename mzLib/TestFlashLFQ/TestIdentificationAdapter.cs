using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;

namespace TestFlashLFQ
{
    internal class TestIdentificationAdapter
    {
        [Test]
        [TestCase(@"FileReadingTests\ExternalFileTypes\FraggerPsm_FragPipev21.1_psm.tsv")]
        public void TestAddProteinGroupInfoCorrect(string path)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, path);
            MsFraggerPsmFile file = new MsFraggerPsmFile(filePath);

            List<Identification> identifications = new List<Identification>();
            identifications = MzLibExtensions.MakeIdentifications(file);

            // list should contain five elements
            Assert.That(identifications.Count, Is.EqualTo(5));
            // one protein associated with given results, list should only contain this one element 
            Assert.That(identifications[0].ProteinGroups.Count, Is.EqualTo(1));
            // two proteins associated with given results, list should contain two elements
            Assert.That(identifications[2].ProteinGroups.Count, Is.EqualTo(2));
            
            Identification identification1= identifications[0];
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

            Dictionary<string, string> allFiles = file.FileNametoFilePath(fullFilePath);

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

            Dictionary<string, string> allFiles = file.FileNametoFilePath(fullFilePath);

            Assert.That(allFiles.TryGetValue(fileName, out var output));
            Assert.AreEqual(output, rawFilePath);
        }
    }
}